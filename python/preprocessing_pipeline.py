# This file contains the preprocessing functions run by run_preprocessing.ipynb

from pathlib import Path
import matplotlib.pyplot as plt
from rt_bids import crop_session, from_session_time
import mne
import numpy as np

# The following 3 functions (_sample_probability, _zscore and jointprob_keep) 
# are python versions of functions from EEGLAB for artifact epoch rejection 
# the code is taken from https://github.com/YoavLev/riseup_preproccesing/blob/yoav/bids/python/Yoav/utils.py 

def _sample_probability(x, bins=1000):
    """For every sample, the probability of its value: a histogram with `bins` bins over the
    range of the data (EEGLAB's realproba)."""
    lowest, highest = x.min(), x.max()
    index = np.floor((x - lowest) / (highest - lowest) * (bins - 1)).astype(int)
    counts = np.bincount(index, minlength=bins)
    return counts[index] / x.size


def _zscore(values):
    """Standardize across epochs (last axis)."""
    return (values - values.mean(axis=-1, keepdims=True)) / values.std(axis=-1, ddof=1, keepdims=True)


def jointprob_keep(epochs_data, threshold):
    """Which epochs to keep, like EEGLAB's pop_jointprob (local and global, threshold in SD).

    epochs_data: array (n_epochs, n_channels, n_times).

    Idea: values that rarely occur in a channel have a low probability. An epoch full of
    improbable values (a large artifact) gets a high "joint probability" score, -sum(log p).
    - local:  the score per channel, standardized across epochs;
    - global: the same score for all channels together.
    An epoch is rejected if its local score in any channel, or its global score, is more
    than `threshold` standard deviations from the mean.
    """
    n_epochs, n_channels, n_times = epochs_data.shape

    local = np.zeros((n_channels, n_epochs))
    for channel in range(n_channels):
        p = _sample_probability(epochs_data[:, channel, :].ravel())
        local[channel] = -np.log(p).reshape(n_epochs, n_times).sum(axis=1)
    local = _zscore(local)

    p_all = _sample_probability(epochs_data.ravel())
    global_ = _zscore(-np.log(p_all).reshape(n_epochs, n_channels * n_times).sum(axis=1))

    reject = (np.abs(local) > threshold).any(axis=0) | (np.abs(global_) > threshold)
    return ~reject


## ICA

def run_ica(raw, sub, ica_output_folder, h_freq = 40, l_freq = 1):
    
    raw_ica = raw.copy().filter(l_freq=l_freq, h_freq=h_freq)

    if raw_ica.info["bads"]:

        raw_ica.interpolate_bads() # CORRMAP requires the ICAs to have the same number of channels. Since some of the raws intiailly contain flagged channels in info["bads"], I have to interpolate them

    raw_ica.set_eeg_reference("average") # rereference, because of interpolation

    raw_ica.resample(250)

    epochs = mne.make_fixed_length_epochs(raw_ica, duration=1, reject_by_annotation=True, preload=True) # making 1s long dummy epochs for ICA
    keep = jointprob_keep(epochs.get_data(picks="eeg"), 3) # using a custom-made function for rejecting dummy epochs that approximates the pop_jointprob in EEGLAB
    epochs_ica = epochs[np.where(keep)[0]]

    ica = mne.preprocessing.ICA(n_components=None, method="infomax", max_iter="auto", rng=97) # when n_components is None, ICA automatically calculates the rank of the data (taking into account the bad channels in info["bads"]) and then decides on the number of components corresponding to it
    ica.fit(epochs_ica, picks="eeg")
    figs = ica.plot_components()

    for i, f in enumerate(figs):
        f.savefig(ica_output_folder / f"{sub}_components_{i+1}.png", dpi=300, bbox_inches='tight')
        plt.close(f)

    ica.save(ica_output_folder / f"{sub}_ica.fif", overwrite=True)
    return ica

# The following functions and code was adapted from Luke Ring's code (https://github.com/NexusDynamic/RiseTogether-analysis/blob/main/explore/s03_ica_explorer.ipynb)

def select_events(session: str, trial_type: str, values, conditions, side: str = "any",group_by: str | None = None) -> pd.DataFrame:
    """The events of one session that the settings ask for, with a "group" column for pooling."""
    ev = events[session]
    chosen = ev[(ev["trial_type"] == trial_type) & ev["value"].isin(values) & ev["condition"].isin(conditions)]
    if side != "any":
        chosen = chosen[chosen["side"] == side]
    column = GROUP_COLUMNS[group_by]
    return chosen.assign(group="all" if column is None else chosen[column].astype(str))

def make_epochs(raw, ev: pd.DataFrame, tmin: float = -0.2, tmax: float = 0.8, baseline: tuple[float | None, float | None] | None = (-0.2, 0.0), picks=("eeg")):
    """Epochs around the events in ev, from a full-session preprocessed `raw` object.. Returns (epochs or None, events in gaps)."""
    k = from_session_time(raw, ev["onset"])  # each event's sample in the recording (MNE's events sample column)
    in_gap = np.isnan(k)                         # the headset skipped samples at that moment
    ev = ev[~in_gap].reset_index(drop=True)
    if ev.empty:
        return None, int(in_gap.sum())
    mne_events = np.column_stack([np.round(k[~in_gap]).astype(int), np.zeros(len(ev), int), np.ones(len(ev), int)])
    epochs = mne.Epochs(raw, mne_events, event_id={"event": 1}, tmin=tmin, tmax=tmax, baseline=baseline,
                        picks=picks, metadata=ev, reject=None, preload=True)
    return epochs, int(in_gap.sum())