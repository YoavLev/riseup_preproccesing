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
BAD_CHANNEL_Z = 5.0  # noisy: robust amplitude more than 5 robust SDs above the median channel (as PREP)
BAD_CHANNEL_FLAT_UV = 1.0  # flat: robust amplitude below 1 µV (1-40 Hz)
BAD_CHANNEL_CORR = 0.4  # uncorrelated: highest correlation with any other channel below 0.4 (as PREP) ...
BAD_CHANNEL_CORR_FRACTION = 0.5  # ... in more than half of the 1 s windows (PREP: 1%; too strict for 27 electrodes)

# Reference of the EEG channels, set in 01_ica after the bad channels are interpolated:
#   None      = keep the recording reference (mastoids), as in the reference pipeline;
#   "average" = common average reference (used, decided 2026-11-11; supervisor's request 2026-11-06).
# The environment variable RT_REFERENCE overrides it.
REFERENCE = os.environ.get("RT_REFERENCE") or "average"

# --- Events (annotations from the dataset's events.tsv, read by rt_bids) ---
CONDITIONS = ["motor_baseline", "individual", "joint"]  # trials to keep (motor baseline kept as a control)
TRIAL_START = "TRIAL_START/"  # e.g. TRIAL_START/joint/r1/t2
TRIAL_END = "TRIAL_END/"
BREAK_MIN_DURATION = 15.0  # gaps between trials longer than this (s) are marked as BAD_break

# Bad stretches (01_ica step 3b, added to the reference pipeline): the EEG (1-40 Hz) is cut into 1 s windows;
# a window is marked BAD_amplitude when the median of its channel SDs is above the threshold.
# Typical clean data: about 10 µV. Meant for stretches where (nearly) all channels are hit at once,
# e.g. headset or cable problems (sub-013, joint block 3); single bad channels are handled in step 4.
BAD_SEGMENT_WINDOW = 1.0  # window length [s]
BAD_SEGMENT_SD_UV = 50.0  # threshold [µV]

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
    # Provisional (session03, 2026-11-11): suggested by scripts/suggest_corrmap_templates.py session03
    # (ICLabel); check the template plots in 02_corrmap, step 2. Replaces the session01 templates
    # (sub-005 4, sub-003 10, sub-004 21), which were fitted with the mastoid reference.
    "eye_blink": ("sub-013", 5),  # ICLabel eye p = 0.99      [ana02: blink_set, blink_comp]
    "lateral_eye": ("sub-016", 9),  # ICLabel eye p = 1.00    [ana02: eyemov_set, eyemov_comp]
    # heartbeat: left out: no component reaches ICLabel heart p >= 0.5 (best: sub-018 8, p = 0.10)
    #   [ana02: heart_set, heart_comp]
    # muscle: not removed (supervisor, 2026-11-11): removing muscle components can also remove brain
    #   signal. (A session03 muscle template, sub-018 5, matched only itself; ICLabel saw 13 more.)
    #   Muscle stays in the data: watch the beta band (central components, tap rate, EMG).
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
# Recorded against the mastoids (confirmed); re-referenced to the common average in 01_ica (REFERENCE above),
# unlike the reference pipeline. M1/M2 themselves are removed from the data (DROP_CHANNELS).

# --- ICA as a spatial filter (04_spatial_filter) ---
# Bands of the envelopes, per area. The game has no sound, so no auditory area.
# Components are chosen by eye (no automatic scoring), see SPATIAL_FILTER_CHOICE.
SPATIAL_FILTER_BANDS = {
    "visual": {"alpha": (8, 12)},
    "motor": {"mu": (8, 12), "beta": (13, 30)},
}
# Chosen after the review with the supervisor: {area: (component, band of the envelope)}, e.g.
# "sub-001": {"visual": (5, "alpha"), "motor": (8, "beta")}. Envelopes are computed only for these.
SPATIAL_FILTER_CHOICE = {
    # "sub-001": {},
}
