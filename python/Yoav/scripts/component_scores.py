"""Numbers per kept session-ICA component, to help choose a visual (alpha) and a motor (mu/beta) component.
Input for component_candidates.py. Prints/saves only aggregate numbers.

Usage (on UCloud, in riseup_preproccesing, after 03_clean and, for a localizer ICA, 04a_localizer_ica):
    uv run python python/Yoav/scripts/component_scores.py session03 [session|visual|motor]
Saves <OUT_ROOT>/component_scores/<session>_<ica>_ica.csv

Per component (time course on the clean data, task time without breaks/bad stretches):
- alpha_db / beta_db: largest peak above the 1/f background (dB) in 8-12 Hz and 13-30 Hz, with its frequency.
  Background = straight line in log-log fitted on 2-7 Hz and 32-40 Hz.
- slope_20_40: log-log slope of the spectrum from 20 to 40 Hz; brain is clearly negative (about -1 to -2),
  muscle is flat or rising (> about -0.5).
- vis_roi / mot_roi: share of the squared map weights in the visual / motor electrodes; top3: strongest electrodes.
- tap_vs_video_alpha/beta_db: power while tapping (motor baseline) minus while watching (video), in dB.
  Motor components: beta (and mu) drop while tapping (negative).
- ICLabel: most likely class and its probability (rough here: ICLabel expects 1-100 Hz, the data are 0.1-40 Hz).
"""
import re
import sys
from pathlib import Path

import mne
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # python/Yoav, for config
import config as cfg  # noqa: E402

VISUAL = ["O1", "O2", "Pz", "P3", "P4", "P7", "P8"]
MOTOR = ["C3", "Cz", "C4", "FC1", "FC2", "CP1", "CP2"]
mne.set_log_level("error")


def window(raw, start_pattern, end_pattern):
    ann = raw.annotations
    offset = raw.first_time if ann.orig_time is not None else 0.0
    t = lambda p: np.array([o - offset for o, d in zip(ann.onset, ann.description) if re.search(p, d)])
    if len(t(start_pattern)) == 0 and "video_first_frame" in start_pattern:   # no photodiode onset (session04)
        start_pattern = r"^MARKER/code \d+ start \(instructions_video\)$"   # 2026-11-24: start marker instead
    start = t(start_pattern).min()
    end = t(end_pattern)[t(end_pattern) > start].min()
    return start, end


def log_psd(raw, reject):
    psd = raw.compute_psd(method="welch", fmin=1, fmax=40, n_fft=1000, n_overlap=0, picks="all",
                          reject_by_annotation=reject)
    return psd.freqs, 10 * np.log10(psd.get_data())  # components x freqs, dB


def band_power(freqs, db, lo, hi):
    return db[:, (freqs >= lo) & (freqs <= hi)].mean(axis=1)


def peaks(freqs, db):
    """Peak above the 1/f line in alpha and beta, and the 20-40 Hz slope, per component."""
    lf = np.log10(freqs)
    fit = ((freqs >= 2) & (freqs <= 7)) | ((freqs >= 32) & (freqs <= 40))
    out = []
    for row in db:
        b = np.polyfit(lf[fit], row[fit], 1)
        resid = row - np.polyval(b, lf)
        res = {}
        for name, lo, hi in [("alpha", 8, 12), ("beta", 13, 30)]:
            m = (freqs >= lo) & (freqs <= hi)
            i = np.argmax(resid[m])
            res[f"{name}_db"] = resid[m][i]
            res[f"{name}_hz"] = freqs[m][i]
        hi = (freqs >= 20) & (freqs <= 40)
        res["slope_20_40"] = np.polyfit(lf[hi], row[hi] / 10, 1)[0]  # log10 power per log10 Hz
        out.append(res)
    return pd.DataFrame(out)


def scores(sub, ica_name="session"):
    out_dir = cfg.OUT_ROOT / sub
    raw = mne.io.read_raw(out_dir / f"{sub}_clean_raw.fif", preload=True)
    if ica_name == "session":
        ica = mne.preprocessing.read_ica(out_dir / f"{sub}_ica.fif")
        removed = pd.read_csv(out_dir / f"{sub}_badcomps.csv")["component"].astype(int).tolist()
    else:  # 2026-11-25: a localizer ICA (04a), fitted on the clean data: every component is kept
        ica = mne.preprocessing.read_ica(out_dir / f"{sub}_{ica_name}_localizer_ica.fif")
        removed = []
    kept = [c for c in range(ica.n_components_) if c not in removed]

    # ICLabel (expects average reference: the data are average-referenced)
    from mne_icalabel import label_components
    labels = label_components(raw, ica, method="iclabel")  # most likely class and its probability

    sources = ica.get_sources(raw)
    freqs, db = log_psd(sources, reject=True)  # task time: breaks and bad stretches left out
    table = peaks(freqs, db)

    # tapping (motor baseline) vs watching (video)
    seg = {}
    for name, (sp, ep) in cfg.LOCALIZERS.items():
        s, e = window(raw, sp, ep)
        f, d = log_psd(ica.get_sources(raw.copy().crop(s, e)), reject=False)  # (sources cannot be cropped in MNE 1.13)
        seg[name] = d
    for band, lo, hi in [("alpha", 8, 12), ("beta", 13, 30)]:
        table[f"tap_vs_video_{band}_db"] = band_power(f, seg["motor"], lo, hi) - band_power(f, seg["visual"], lo, hi)

    maps = ica.get_components() ** 2
    maps = maps / maps.sum(axis=0)
    names = ica.info["ch_names"]
    table["vis_roi"] = maps[[names.index(c) for c in VISUAL]].sum(axis=0)
    table["mot_roi"] = maps[[names.index(c) for c in MOTOR]].sum(axis=0)
    table["top3"] = [" ".join(names[i] for i in np.argsort(maps[:, k])[::-1][:3]) for k in range(maps.shape[1])]

    table["iclabel"] = labels["labels"]
    table["iclabel_p"] = labels["y_pred_proba"]
    table.insert(0, "component", range(ica.n_components_))
    table.insert(0, "subject", sub)
    return table[table["component"].isin(kept)]


if __name__ == "__main__":
    session = sys.argv[1]
    ica_name = sys.argv[2] if len(sys.argv) > 2 else "session"   # session | visual | motor (localizer ICA)
    result = pd.concat([scores(s, ica_name) for s in cfg.session_subjects(session)])
    out = cfg.OUT_ROOT / "component_scores" / f"{session}_{ica_name}_ica.csv"
    out.parent.mkdir(exist_ok=True)
    result.round(3).to_csv(out, index=False)
    print("saved", out)
