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

# --- Filters ---
# All filtering uses MNE's filter() (default FIR design: the filter length follows from the cut-off).
# Not the reference's fixed filter orders (ana01/ana03: LP_ord 100, HP_ord 500; EEGLAB pop_firws, Hann
# window): the recordings drift slowly by several mV, and a 1 Hz high-pass only 2 s long (order 500 at
# 250 Hz) lets ~1/80 of that drift through (tens of µV, more than the EEG); the 0.1 Hz high-pass of
# order 500 hardly filters at all. MNE's filters are long enough. See ica_reference_questions.md #9.

# --- ICA training copy (ana01_ICA.m) ---
ICA_LP = 40  # low-pass cut-off [Hz], only for ICA                      [ana01: LP]
ICA_SRATE = 250  # downsample the copy for ICA [Hz]                      [ana01: SRATE]
ICA_HP = 1  # high-pass cut-off [Hz], only for ICA                       [ana01: HP]
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
    # One set per session (02_corrmap matches each session to its own templates).
    # 2026-11-24: ICAs refitted with the dataset's bad channels (channels.tsv, Luke's update): component numbers changed.
    # Chosen as before: suggested by scripts/suggest_corrmap_templates.py <session> (ICLabel), checked on the template
    # plots in 02_corrmap, step 2. (Before 2026-11-24, session03: eye_blink sub-013 1, lateral_eye sub-014 1.)
    # heartbeat: [ana02: heart_set, heart_comp]; left out so far (no clear heart component).
    # muscle: not removed (supervisor, 2026-11-11): removing muscle components can also remove brain
    #   signal. Muscle stays in the data: watch the beta band (central components, tap rate, EMG).
    "session03": {  # 2026-11-24: ICLabel suggestions, maps checked by eye
        "eye_blink": ("sub-013", 0),  # eye p = 1.00, symmetric frontal map   [ana02: blink_set, blink_comp]
        "lateral_eye": ("sub-014", 1),  # eye p = 1.00, F7/F8 opposite signs    [ana02: eyemov_set, eyemov_comp]
    },
    "session04": {  # 2026-11-24: ICLabel suggestions, maps checked by eye
        "eye_blink": ("sub-024", 0),  # eye p = 1.00, symmetric frontal map
        "lateral_eye": ("sub-019", 1),  # eye p = 1.00, F7/F8 opposite signs
    },
    "session02": {  # 2026-11-25: ICLabel suggestions, maps checked by eye
        "eye_blink": ("sub-010", 0),  # eye p = 1.00, symmetric frontal map
        "lateral_eye": ("sub-009", 2),  # eye p = 1.00, F7/F8 opposite signs
    },
    "session05": {  # 2026-11-25: ICLabel suggestions, maps checked by eye
        "eye_blink": ("sub-027", 0),  # eye p = 1.00, symmetric frontal map
        "lateral_eye": ("sub-026", 1),  # eye p = 1.00, F7/F8 opposite signs
    },
    "session06": {  # 2026-11-25: ICLabel suggestions, maps checked by eye
        "eye_blink": ("sub-032", 0),  # eye p = 0.99, symmetric frontal map
        "lateral_eye": ("sub-031", 2),  # eye p = 1.00, F7/F8 opposite signs
    },
}
CORRMAP_THRESHOLD = 0.8  # correlation threshold                         [ana02: 'th', '0.8']
# Extra components to remove per participant after visual review (beyond what CORRMAP found).
# Rule (as on 2026-11-12): ICLabel eye p >= 0.8 and the map fits neither template; checked by eye.
# (Before 2026-11-24, old ICAs: sub-013 [7], sub-015 [8], sub-018 [13]; numbers no longer valid.)
EXTRA_EXCLUDE = {
    # 2026-11-24 (new ICAs): ICLabel eye p >= 0.8, missed by CORRMAP; frontal-pole maps checked by eye
    "sub-013": [8],   # eye p = 0.90
    "sub-018": [13],  # eye p = 0.94
    # 2026-11-25, sessions 02/05/06 (same rule; all frontal-pole maps, checked by eye)
    "sub-008": [23],  # eye p = 0.94
    "sub-011": [9],   # eye p = 0.81 (borderline; frontal band Fp1-Fp2-F7-F8)
    "sub-012": [15],  # eye p = 0.99
    "sub-025": [9],   # eye p = 0.97
    "sub-027": [15],  # eye p = 0.98
    "sub-028": [19],  # eye p = 0.93
    "sub-029": [19],  # eye p = 0.91
    "sub-032": [20],  # eye p = 0.94
    "sub-033": [13],  # eye p = 0.98
    "sub-035": [20],  # eye p = 0.95
    "sub-036": [17],  # eye p = 0.92
}
# Components CORRMAP found that are kept anyway, after visual review (false matches)
# (Before 2026-11-24: sub-017 [16], matched the blink template but not a blink.)
EXTRA_KEEP = {
    "sub-017": [15],  # 2026-11-24: matched the blink template, but ICLabel eye p = 0.05 (other 0.81); as on 2026-11-12
    # 2026-11-25, sessions 02/05/06: matched the blink template, but ICLabel eye p < 0.1 and a broad frontal-central
    # dipole map rather than a frontal-pole one (same rule as sub-017)
    "sub-007": [6],   # brain 1.00
    "sub-008": [19],  # other; eye p = 0.07
    "sub-009": [9],   # brain 1.00
    "sub-012": [2, 13],  # brain 0.99 / 0.97; ICA002 is the most blink-like of these (Yoav to check)
    "sub-026": [8],   # brain 0.99 (frontal-midline maximum)
    "sub-029": [11],  # brain 1.00
}

# --- Analysis data (ana03_preprocessing.m) ---
LP = 40  # low-pass cut-off [Hz]                                         [ana03: LP_erp]
HP = 0.1  # high-pass cut-off [Hz]                                       [ana03: HP_erp]
# Recorded against the mastoids (confirmed); re-referenced to the common average in 01_ica (REFERENCE above),
# unlike the reference pipeline. M1/M2 themselves are removed from the data (DROP_CHANNELS).

# --- Visual localizer ICA (04a_localizer_ica; supervisor's request, 2026-11-18) ---
# Only occipital alpha is analysed (decided 2026-11-28; no motor components). A separate ICA on the visual
# localizer, the instruction video (watched without playing), fitted on the clean data (03_clean, eye components
# already removed), with the same settings as 01_ica (ICA_HP, ICA_LP, ICA_SRATE, ICA_EPOCH_LENGTH, ICA_PRUNE,
# RANDOM_STATE). One occipital alpha component of it is then used as a spatial filter on the game trials: the filter
# is learned from data the analysis does not test.
# Window: (regular expression of the start annotation, of the end annotation), searched in each recording:
# from the video's first frame on this participant's screen (photodiode) to the end of the phase.
LOCALIZER = (r"^DISPLAY_ONSET/video_first_frame$", r"^MARKER/code \d+ stop \(instructions_video\)$")
# Start used when a recording has no photodiode onset (session04 has no DISPLAY_ONSET/video_first_frame; 2026-11-25).
# The phase marker is sent when the video is started, so the window may begin slightly before the first frame.
LOCALIZER_START_FALLBACK = r"^MARKER/code \d+ start \(instructions_video\)$"
LOCALIZER_MARGIN = 10.0  # s of data kept on each side while filtering (no filter edge effects), then cut off

# --- Choosing the occipital alpha component (scripts/component_scores.py, component_candidates.py) ---
# Shortlist rule; the final choice is made by eye on the candidate PNGs. A component is a candidate if
#   - its strongest electrode is occipital (OCCIPITAL; the cap has no Oz),
#   - at least POSTERIOR_MIN_SHARE of its map's squared weights lie in the posterior electrodes (POSTERIOR),
#   - its spectrum is brain-like: 20-40 Hz log-log slope below MAX_SLOPE_20_40 (muscle is flat or rising),
#   - it has an alpha peak at least MIN_ALPHA_DB above the 1/f background;
# candidates are ranked by the alpha peak. A participant without a candidate is excluded.
OCCIPITAL = ["O1", "O2"]
POSTERIOR = ["O1", "O2", "P7", "P8", "P3", "P4", "Pz"]
POSTERIOR_MIN_SHARE = 0.35
MAX_SLOPE_20_40 = -0.5
MIN_ALPHA_DB = 1.0

# --- ICA as a spatial filter (04_spatial_filter) ---
ALPHA_BAND = (8, 12)  # band of the envelope [Hz]
# Chosen by eye from the candidate PNGs: {participant: component of the visual localizer ICA}, e.g. "sub-015": 6.
# Envelopes are computed only for these; participants missing here are excluded (say why in a comment).
SPATIAL_FILTER_CHOICE = {
}
# Width of the band-pass slopes [Hz], both sides. MNE's default is wider (2 Hz at 8 Hz), which lets theta and beta
# leak into the alpha envelope; with 1 Hz the -6 dB points are 7.5 and 12.5 Hz. The filter is then 3.3 s long, no
# problem on the continuous data.
SPATIAL_FILTER_TRANSITION = 1.0
# Sampling rate of the saved envelopes [Hz]. cRQA compares every time point with every other one, so its cost grows
# with the square of the length: a 300 s trial has 150 000 points at 500 Hz but 7 500 at 25 Hz.
# The envelopes change slowly, so little is lost (resampling low-passes at half this rate, 12.5 Hz).
ENVELOPE_SFREQ = 25

# --- cRQA (05_windows.ipynb cuts the windows, 06_crqa.R runs the crqa R package) ---
# Windows and results go to their own folder, next to preproc: derivatives/crqa/<session>/
# (the cRQA settings are at the top of 06_crqa.R)
CRQA_ROOT = OUT_ROOT.parent / "crqa"
# Smoothing of the envelopes before cRQA (05_windows): low-pass at this frequency [Hz], or None (no smoothing).
# Zero-phase FIR filter with a Hamming window (MNE filter(), firwin design; slope half the cut-off).
# A narrow band (alpha 8-12 Hz) gives an envelope with fast bumps of ~0.3-1 s (it "beats" at about the band's
# width) that are unlikely to be shared between two people; a 1 Hz low-pass keeps only the slow ups and downs
# (see 04, section 3b, panel 4). Applied to the whole 25 Hz envelope before the windows are cut (no edge effects).
# With smoothing, the windows go to their own folder, derivatives/crqa/<session>_smooth<X>hz/, so both can be
# compared. Not decided yet (2026-11-24: tried for the demo). The environment variable RT_SMOOTH_HZ overrides it
# ("0" = no smoothing).
# Decided 2026-11-25 (after the comparison of 0.5, 1 and 2 Hz): 1 Hz.
CRQA_SMOOTH_HZ = 1.0
# Sampling rate of the windows after the smoothing [Hz], or None (keep ENVELOPE_SFREQ, 25 Hz). After a low-pass at
# f Hz (slope f/2: nothing left above 1.5 f), 5 f samples per second keep everything (Nyquist 2.5 f) and drop the
# redundant samples: at 25 Hz a 1 Hz envelope barely changes from one sample to the next, so almost every recurrent
# point lies on a diagonal line and %DET is ~98 % for any pair (demo, 2026-11-24). Rule used: 5 x CRQA_SMOOTH_HZ.
# Output folder then: derivatives/crqa/<session>_smooth<X>hz_<Y>hz/. RT_CRQA_SFREQ overrides it.
# Decided 2026-11-25: 5 Hz (with CRQA_SMOOTH_HZ = 1).
CRQA_SFREQ = 5.0
