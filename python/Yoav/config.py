"""Settings for the three preprocessing notebooks (01_ica, 02_corrmap, 03_clean).

The ICA and filtering constants follow the reference pipeline the lab uses:
Stropahl, Bauer, Debener & Bleichner (2018), "Source-modelling auditory processes of EEG data
using EEGLAB and Brainstorm", supplementary material (notes/papers/eeg_data_sheet.PDF).
The name of each constant in the reference scripts is given in brackets, e.g. [ana01: LP].
Deviations and open questions are explained in ica_reference_questions.md.
"""

import os
import sys
from pathlib import Path

# --- Paths ---
# The BIDS dataset (read-only, on UCloud only: the data never leave UCloud).
# RT_BIDS_ROOT overrides it, e.g. for a test copy.
# The dataset is mounted either as its parent folder or as the dataset folder itself, depending on the job
_BIDS_MOUNTS = [Path("/work/RiseTogetherEEG/VideoGameParadigm-EEG"), Path("/work/VideoGameParadigm-EEG")]
BIDS_ROOT = Path(os.environ.get("RT_BIDS_ROOT") or next((p for p in _BIDS_MOUNTS if p.exists()), _BIDS_MOUNTS[0]))
# All outputs go to my own folder (the dataset is read-only); one folder per participant.
OUT_ROOT = Path(os.environ.get("RT_OUT_ROOT", "/work/YoavLevy#3923/derivatives/preproc"))
LOG_DIR = OUT_ROOT.parent.parent / "logs"  # executed notebooks and script output

# The reader shipped with the dataset (code/rt_bids.py); on my Mac (no data) the copy in
# RiseTogether-analysis/bids_code is used, so the notebooks can at least be imported there.
RT_BIDS_CODE = next(p for p in [BIDS_ROOT / "code", Path("~/Bsc/RiseTogether-analysis/bids_code").expanduser()]
                    if (p / "rt_bids.py").exists())
sys.path.insert(0, str(RT_BIDS_CODE))
TASK = "risetogether"


def out_dir(subject):
    """Output folder of one participant, e.g. .../derivatives/preproc/sub-001/ (created if needed)."""
    path = OUT_ROOT / subject
    path.mkdir(parents=True, exist_ok=True)
    return path


def bdf_path(subject):
    """The participant's EEG recording in the dataset."""
    return BIDS_ROOT / subject / "eeg" / f"{subject}_task-{TASK}_eeg.bdf"


def session_subjects(group):
    """The participants (sub-NNN) of one session (group in participants.tsv, e.g. "session01")."""
    import rt_bids
    return rt_bids.session_participants(BIDS_ROOT, group)["participant_id"].tolist()


def session_of(subject):
    """The session (group) a participant belongs to."""
    import pandas as pd
    people = pd.read_csv(BIDS_ROOT / "participants.tsv", sep="\t", na_values="n/a")
    return people.set_index("participant_id").loc[subject, "group"]


# --- Channels ---
MONTAGE = "colin27_1020"  # electrode positions (the new name of "standard_1020" in MNE 1.13)
EOG_CHANNELS = ["EOG2", "EOG3", "EOG4"]  # eye electrodes (there is no EOG1)
DROP_CHANNELS = ["M1", "M2"]  # mastoids: removed from the data (supervisor's request)

# Bad channels per participant, filled in by hand after inspecting the data (01_ica, step 4).
# Channels marked bad in the dataset (channels.tsv status = bad: at the ADC limit > 50% of the
# time) are added automatically; list only the extra ones here. Keys: participant_id (sub-NNN).
BAD_CHANNELS = {
    # "sub-001": [],
}
# Clearly bad channels are also found automatically (utils.find_bad_channels, criteria of the PREP
# pipeline); the thresholds:
BAD_CHANNEL_Z = 5.0  # noisy: SD more than 5 robust SDs above the median channel
BAD_CHANNEL_FLAT_UV = 1.0  # flat: SD below 1 µV (1-40 Hz)
BAD_CHANNEL_CORR = 0.4  # uncorrelated: highest correlation with any other channel below 0.4 ...
BAD_CHANNEL_CORR_FRACTION = 0.01  # ... in more than 1% of the 1 s windows

# --- Events (annotations from the dataset's events.tsv, read by rt_bids) ---
CONDITIONS = ["motor_baseline", "individual", "joint"]  # trials to keep (motor baseline kept as a control)
TRIAL_START = "TRIAL_START/"  # e.g. TRIAL_START/joint/r1/t2
TRIAL_END = "TRIAL_END/"
BREAK_MIN_DURATION = 15.0  # gaps between trials longer than this (s) are marked as BAD_break

# --- ICA training copy (ana01_ICA.m) ---
ICA_LP = 40  # low-pass cut-off [Hz], only for ICA                      [ana01: LP]
ICA_LP_ORDER = 100  # low-pass filter order                              [ana01: LP_ord]
ICA_SRATE = 250  # downsample the copy for ICA [Hz]                      [ana01: SRATE]
ICA_HP = 1  # high-pass cut-off [Hz], only for ICA                       [ana01: HP]
ICA_HP_ORDER = 500  # high-pass filter order                             [ana01: HP_ord]
ICA_EPOCH_LENGTH = 1.0  # dummy epochs of 1 s                            [ana01: eeg_regepochs 'recurrence', 1]
ICA_PRUNE = 3  # reject dummy epochs with joint probability > 3 SD       [ana01: PRUNE]
ICA_PCA = True  # reduce dimensions with PCA before ICA                  [ana01: PCA]
ICA_PCADIMS = 50  # PCA dimensions (capped at the data rank, see questions) [ana01: PCADIMS]
RANDOM_STATE = 97  # fixed seed so ICA gives the same result every run

# --- CORRMAP (ana02_corrmap.m) ---
# Templates: (participant, component) for each artifact type, chosen by hand from the
# component plots of all participants (01_ica saves them as PNG).
# Note: MNE counts components from 0, EEGLAB from 1 (EEGLAB component 1 = MNE component 0).
CORRMAP_TEMPLATES = {
    "eye_blink": ("sub-001", 0),  # TODO: choose after looking at the component plots   [ana02: blink_set, blink_comp]
    "heartbeat": ("sub-001", 1),  # TODO                                                  [ana02: heart_set, heart_comp]
    "lateral_eye": ("sub-001", 2),  # TODO                                                [ana02: eyemov_set, eyemov_comp]
    # Added to the reference (supervisor's decision): tapping produces muscle artifacts in the beta band.
    # Muscle components are focal and often one-sided; if one template misses the other side,
    # add a second one (e.g. "muscle_left" and "muscle_right").
    "muscle": ("sub-001", 3),  # TODO
}
CORRMAP_THRESHOLD = 0.8  # correlation threshold                         [ana02: 'th', '0.8']
# Extra components to remove per participant after visual review (beyond what CORRMAP found)
EXTRA_EXCLUDE = {
    # "sub-001": [],
}

# --- Analysis data (ana03_preprocessing.m) ---
LP = 40  # low-pass cut-off [Hz]                                         [ana03: LP_erp]
LP_ORDER = 100  # low-pass filter order                                  [ana03: LP_ord]
HP = 0.1  # high-pass cut-off [Hz]                                       [ana03: HP_erp]
HP_ORDER = 500  # high-pass filter order                                 [ana03: HP_ord]
# No re-referencing: the data stay referenced to the mastoids, the recording reference (confirmed),
# as in the reference pipeline. M1/M2 themselves are removed from the data (DROP_CHANNELS).

# --- ICA as a spatial filter (04_spatial_filter) ---
# A component is a candidate for an area when its map is mostly over the area's electrodes and its
# spectrum has a peak in the area's band(s). The game has no sound, so no auditory area.
SPATIAL_FILTER = {
    "visual": {"region": ["O1", "O2", "Pz", "P3", "P4", "P7", "P8"], "bands": {"alpha": (8, 12)}},
    "motor": {"region": ["C3", "Cz", "C4", "FC1", "FC2", "CP1", "CP2"], "bands": {"mu": (8, 12), "beta": (13, 30)}},
}
MOTOR_LEFT, MOTOR_RIGHT = ["C3", "FC1", "CP1"], ["C4", "FC2", "CP2"]  # to note which side a motor component is on
SPATIAL_FILTER_BACKGROUND = [(2, 7), (31, 40)]  # frequencies (Hz) used to fit the 1/f background
SPATIAL_FILTER_MIN_LOCATION = 0.4  # share of the component's map weight that must be in the region
SPATIAL_FILTER_N_CANDIDATES = 3  # candidates per area saved as PNG for the review
# Chosen after the review with the supervisor: {area: (component, band of the envelope)}, e.g.
# "sub-001": {"visual": (5, "alpha"), "motor": (8, "beta")}. Envelopes are computed only for these.
SPATIAL_FILTER_CHOICE = {
    # "sub-001": {},
}
