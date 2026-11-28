"""Numbers per component of the visual localizer ICA (04a_localizer_ica), to help choose the occipital alpha component.
Input for component_candidates.py. Prints/saves only aggregate numbers.

Usage (on UCloud, in riseup_preproccesing, after 04a_localizer_ica for every participant of the session):
    uv run python python/Yoav/scripts/component_scores.py session03
Saves <OUT_ROOT>/component_scores/<session>_visual_ica.csv

Per component (time course on the clean data, task time without breaks/bad stretches):
- alpha_db / alpha_hz: largest peak above the 1/f background (dB) in ALPHA_BAND (8-12 Hz), and its frequency.
  Background = straight line in log-log fitted on 2-7 Hz and 32-40 Hz.
- slope_20_40: log-log slope of the spectrum from 20 to 40 Hz; brain is clearly negative (about -1 to -2),
  muscle is flat or rising (> about -0.5).
- posterior_share: share of the squared map weights in the posterior electrodes (POSTERIOR in config.py);
  top3: the three strongest electrodes.
- ICLabel: most likely class and its probability (rough here: ICLabel expects 1-100 Hz, the data are 0.1-40 Hz).
"""
import sys
from pathlib import Path

import mne
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # python/Yoav, for config
import config as cfg  # noqa: E402

mne.set_log_level("error")


def one_over_f_line(freqs, db):
    """The 1/f background of a spectrum in dB: straight line in log-log fitted on 2-7 Hz and 32-40 Hz."""
    fit = ((freqs >= 2) & (freqs <= 7)) | ((freqs >= 32) & (freqs <= 40))
    return np.polyval(np.polyfit(np.log10(freqs[fit]), db[fit], 1), np.log10(freqs))


def scores(sub):
    out_dir = cfg.out_dir(sub)
    raw = mne.io.read_raw(out_dir / f"{sub}_clean_raw.fif", preload=True)
    ica = mne.preprocessing.read_ica(out_dir / f"{sub}_visual_localizer_ica.fif")

    # ICLabel (expects average reference: the data are average-referenced)
    from mne_icalabel import label_components
    labels = label_components(raw, ica, method="iclabel")

    # spectrum of every component, task time: breaks and bad stretches left out
    psd = ica.get_sources(raw).compute_psd(method="welch", fmin=1, fmax=40, n_fft=1000, n_overlap=0, picks="all",
                                           reject_by_annotation=True)
    freqs, db = psd.freqs, 10 * np.log10(psd.get_data())   # components x freqs, dB
    in_band = (freqs >= cfg.ALPHA_BAND[0]) & (freqs <= cfg.ALPHA_BAND[1])
    high = (freqs >= 20) & (freqs <= 40)
    rows = []
    for row in db:
        above = (row - one_over_f_line(freqs, row))[in_band]
        i = np.argmax(above)
        rows.append({"alpha_db": above[i], "alpha_hz": freqs[in_band][i],
                     "slope_20_40": np.polyfit(np.log10(freqs[high]), row[high] / 10, 1)[0]})  # log10 power per log10 Hz
    table = pd.DataFrame(rows)

    maps = ica.get_components() ** 2
    maps = maps / maps.sum(axis=0)
    names = ica.info["ch_names"]
    table["posterior_share"] = maps[[names.index(c) for c in cfg.POSTERIOR]].sum(axis=0)
    table["top3"] = [" ".join(names[i] for i in np.argsort(maps[:, k])[::-1][:3]) for k in range(maps.shape[1])]
    table["iclabel"] = labels["labels"]
    table["iclabel_p"] = labels["y_pred_proba"]
    table.insert(0, "component", range(ica.n_components_))
    table.insert(0, "subject", sub)
    return table


if __name__ == "__main__":
    session = sys.argv[1]
    result = pd.concat([scores(s) for s in cfg.session_subjects(session)])
    out = cfg.OUT_ROOT / "component_scores" / f"{session}_visual_ica.csv"
    out.parent.mkdir(exist_ok=True)
    result.round(3).to_csv(out, index=False)
    print("saved", out)
