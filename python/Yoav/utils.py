"""Python versions of the two EEGLAB functions used in the reference pipeline
(Stropahl et al., 2018, supplementary material): pop_firws (filtering) and pop_jointprob
(rejecting unusual epochs before ICA).
"""

import numpy as np
from scipy.signal import firwin, oaconvolve


def firws_filter(raw, cutoff, kind, order):
    """Filter like EEGLAB's pop_firws: windowed-sinc FIR filter with a Hann window.

    cutoff: the cut-off frequency in Hz (the -6 dB point).
    kind:   "lowpass" or "highpass".
    order:  filter order (even number); the filter has order + 1 coefficients.

    The filter is applied once and is zero-phase (its delay is compensated), so events stay
    at the same time. EEG and EOG channels are filtered in place.
    """
    if order % 2:
        raise ValueError("the filter order must be even")
    coefficients = firwin(order + 1, cutoff, window="hann", pass_zero=(kind == "lowpass"),
                          fs=raw.info["sfreq"])
    half = order // 2

    def apply(signal):
        # Repeat the first/last sample at the edges, so the edges are not pulled towards zero
        padded = np.pad(signal, half, mode="edge")
        return oaconvolve(padded, coefficients, mode="valid")

    raw.apply_function(apply, picks=["eeg", "eog"], channel_wise=True)
    return raw


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


def find_bad_channels(raw, z_max=5.0, corr_min=0.4, corr_fraction=0.5, flat_uv=1.0, eye_channels=("Fp1", "Fp2"),
                      window=1.0):
    """Clearly bad EEG channels, with the two main criteria of the PREP pipeline (Bigdely-Shamlo et al.,
    2015), on a 1-40 Hz copy of the task data (breaks left out).

    - noisy or flat: the channel's robust amplitude (0.7413 x interquartile range, as PREP; rare
      large events such as blinks hardly change it) is more than `z_max` robust SDs above the median
      channel, or below `flat_uv` microvolts. The eye-near channels (`eye_channels`) are not checked
      for noise: their large amplitude is eye activity, which ICA removes (CORRMAP's eye templates
      need these channels);
    - not correlated: in each 1 s window, the channel's highest absolute correlation with any other
      channel; bad if it is below `corr_min` in more than `corr_fraction` of the windows. PREP uses
      1% of the windows; with only 27 electrodes the temporal and frontal channels have no close
      neighbour and often fall below, so here a channel must be uncorrelated most of the time.

    Returns a table with one row per channel and the reason it is bad ("" if it is good).
    """
    import mne
    import pandas as pd

    eeg = raw.copy().pick(mne.pick_types(raw.info, eeg=True, exclude=[]))
    eeg.filter(1, 40, verbose="error")
    x = eeg.get_data(reject_by_annotation="omit")

    q75, q25 = np.percentile(x, [75, 25], axis=1)
    amplitude = 0.7413 * (q75 - q25)
    z = (amplitude - np.median(amplitude)) / (1.4826 * np.median(np.abs(amplitude - np.median(amplitude))))

    n = int(window * raw.info["sfreq"])
    w = x[:, : x.shape[1] // n * n].reshape(len(x), -1, n).transpose(1, 0, 2)  # (windows, channels, samples)
    w = w - w.mean(axis=2, keepdims=True)
    w = w / (w.std(axis=2, keepdims=True) + 1e-20)
    corr = np.abs(np.einsum("wcs,wds->wcd", w, w) / n)
    corr[:, np.arange(len(x)), np.arange(len(x))] = 0  # ignore each channel's correlation with itself
    low_corr = (corr.max(axis=2) < corr_min).mean(axis=0)

    table = pd.DataFrame({"channel": eeg.ch_names, "amplitude_uv": amplitude * 1e6, "amplitude_z": z,
                          "low_corr_fraction": low_corr})
    reasons = []
    for row in table.itertuples():
        reason = []
        if row.amplitude_z > z_max and row.channel not in eye_channels:
            reason.append("noisy")
        if row.amplitude_uv < flat_uv:
            reason.append("flat")
        if row.low_corr_fraction > corr_fraction:
            reason.append("uncorrelated")
        reasons.append(" ".join(reason))
    table["reason"] = reasons
    return table.round(3)
