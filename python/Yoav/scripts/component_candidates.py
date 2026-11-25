"""Shortlist of visual (alpha) and motor (beta) components per participant, for choosing SPATIAL_FILTER_CHOICE by eye.

Input: the numbers of /work/YoavLevy#3923/tests/component_scores.py <session> (per kept session-ICA component:
alpha/beta peak above the 1/f background, 20-40 Hz slope, share of the map in the visual/motor electrodes, strongest
electrodes). The rule written down from how the session03 trial components were chosen (2026-11-19; it gives the
same component for 3 of the 4, see below), and from CLAUDE.md ("one alpha component (parietal/occipital, peak above
1/f) and one beta component (central)"):

  visual: strongest electrode in the visual ROI (O1 O2 Pz P3 P4 P7 P8), at least 35 % of the map's squared weights
          in that ROI, brain-like spectrum (20-40 Hz slope below -0.5; muscle is flat), ranked by the alpha peak (dB);
  motor:  the same with the motor ROI (C3 Cz C4 FC1 FC2 CP1 CP2), ranked by the beta peak (dB).

(2026-11-19 choices on the old ICAs: sub-015 motor 13, sub-017 visual 7 and motor 19 are the rule's first; sub-015
visual 4 is not: its 20-40 Hz slope was -0.12, flat like muscle.)

The top 3 of each area are drawn (scalp map + spectrum with the 1/f line and the band) so the final choice is made by
eye, as before. Saves <OUT_ROOT>/spatial_filter_review/<session>/candidates_<sub>.png and candidates_<session>.csv.

Localizer ICA (2026-11-25): with a second argument (visual or motor) the candidates come from that localizer ICA
(04a_localizer_ica; scores from component_scores.py <session> <localizer>), only for that area. A fourth column then
shows the current choice from the session ICA (SPATIAL_FILTER_CHOICE) for comparison, and each candidate gets
map_r (|correlation| of its scalp map with that choice's map) and env_r (correlation of their 8-12 / 13-30 Hz
Hilbert envelopes on the clean data, bad stretches left out): values near 1 mean the same source.
Saves candidates_<sub>_<localizer>.png and candidates_<session>_<localizer>.csv.

Usage (on UCloud, after 03_clean and component_scores.py):
    uv run python python/Yoav/scripts/component_candidates.py session03            # session ICA, both areas
    uv run python python/Yoav/scripts/component_candidates.py session03 visual     # visual localizer ICA
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as cfg  # noqa: E402

mne.set_log_level("error")
SESSION = sys.argv[1]
ICA_NAME = sys.argv[2] if len(sys.argv) > 2 else "session"   # session | visual | motor
SCORE_DIR = Path("/work/YoavLevy#3923/tests/component_scores")
SCORES = SCORE_DIR / f"{SESSION}_{ICA_NAME}_ica.csv"
OUT = cfg.OUT_ROOT / "spatial_filter_review" / SESSION
SUFFIX = "" if ICA_NAME == "session" else f"_{ICA_NAME}"
ROI = {"visual": ["O1", "O2", "Pz", "P3", "P4", "P7", "P8"], "motor": ["C3", "Cz", "C4", "FC1", "FC2", "CP1", "CP2"]}
if ICA_NAME != "session":
    ROI = {ICA_NAME: ROI[ICA_NAME]}   # a localizer ICA is searched only for its own area
RANK_BY = {"visual": "alpha_db", "motor": "beta_db"}
BAND = {"visual": (8, 12), "motor": (13, 30)}
MIN_ROI_SHARE = 0.35
MAX_SLOPE = -0.5
N_SHOW = 3

scores = pd.read_csv(SCORES)
scores["top1"] = scores["top3"].str.split().str[0]
rows = []
for area, roi in ROI.items():
    share = "vis_roi" if area == "visual" else "mot_roi"
    ok = scores["top1"].isin(roi) & (scores[share] >= MIN_ROI_SHARE) & (scores["slope_20_40"] < MAX_SLOPE)
    ranked = scores[ok].sort_values(["subject", RANK_BY[area]], ascending=[True, False])
    ranked = ranked.assign(area=area, rank=ranked.groupby("subject").cumcount() + 1)
    rows.append(ranked[ranked["rank"] <= N_SHOW])
cand = pd.concat(rows)[["subject", "area", "rank", "component", "alpha_db", "alpha_hz", "beta_db", "beta_hz",
                        "slope_20_40", "vis_roi", "mot_roi", "top3", "tap_vs_video_beta_db", "iclabel", "iclabel_p"]]


def ica_path(sub, name):
    return cfg.out_dir(sub) / (f"{sub}_ica.fif" if name == "session" else f"{sub}_{name}_localizer_ica.fif")


def band_envelope(ica, raw, comp, band):
    """Hilbert envelope of one component in the band, bad stretches left out (as 04, without downsampling)."""
    src = ica.get_sources(raw).pick([f"ICA{comp:03d}"])
    src.filter(*band, picks="all", l_trans_bandwidth=cfg.SPATIAL_FILTER_TRANSITION,
               h_trans_bandwidth=cfg.SPATIAL_FILTER_TRANSITION).apply_hilbert(envelope=True, picks="all")
    return src.get_data(reject_by_annotation="omit")[0]


reference = {}   # subject -> (session ICA, component) of the current choice, for a localizer ICA
if ICA_NAME != "session":
    session_scores = pd.read_csv(SCORE_DIR / f"{SESSION}_session_ica.csv")
    for sub in cfg.session_subjects(SESSION):
        if ICA_NAME in cfg.SPATIAL_FILTER_CHOICE.get(sub, {}):
            ica_name, comp, _ = cfg.SPATIAL_FILTER_CHOICE[sub][ICA_NAME]
            reference[sub] = (ica_name, comp, session_scores[(session_scores["subject"] == sub)
                                                             & (session_scores["component"] == comp)].iloc[0])
    cand["map_r"] = np.nan
    cand["env_r"] = np.nan

for sub in cfg.session_subjects(SESSION):
    c = cand[cand["subject"] == sub]
    if c.empty:
        continue
    raw = mne.io.read_raw(cfg.out_dir(sub) / f"{sub}_clean_raw.fif", preload=True)
    ica = mne.preprocessing.read_ica(ica_path(sub, ICA_NAME))
    ref = reference.get(sub)
    n_col = N_SHOW + (ref is not None)
    fig, axes = plt.subplots(len(ROI), 2 * n_col, figsize=(4.2 * n_col, 3.2 * len(ROI) + 0.4),
                             gridspec_kw={"width_ratios": [1, 1.4] * n_col}, squeeze=False)
    if ref is not None:
        ref_ica = mne.preprocessing.read_ica(ica_path(sub, ref[0]))
        ref_map = ref_ica.get_components()[:, ref[1]]
        ref_env = band_envelope(ref_ica, raw, ref[1], BAND[ICA_NAME])

    def draw(ax_map, ax_psd, this_ica, comp, r, area, title):
        mne.viz.plot_topomap(this_ica.get_components()[:, comp], this_ica.info, axes=ax_map, show=False, contours=0)
        ax_map.set_title(title, fontsize=9)
        src = this_ica.get_sources(raw).pick([f"ICA{comp:03d}"])
        psd = src.compute_psd(fmin=1, fmax=40, picks="all", n_fft=1000, reject_by_annotation=True)
        f, p = psd.freqs, 10 * np.log10(psd.get_data()[0])
        fit = ((f >= 2) & (f <= 7)) | ((f >= 32) & (f <= 40))     # 1/f line, as in component_scores.py
        b = np.polyfit(np.log10(f[fit]), p[fit], 1)
        ax_psd.plot(f, p, color="k", lw=1)
        ax_psd.plot(f, np.polyval(b, np.log10(f)), color="0.6", lw=1, ls="--")
        ax_psd.axvspan(*BAND[area], color="C0", alpha=0.15)
        peak = RANK_BY[area]
        ax_psd.set_title(f"{peak.split('_')[0]} peak {r[peak]:.1f} dB at {r[peak.replace('db', 'hz')]:g} Hz; "
                         f"slope {r['slope_20_40']:.1f}", fontsize=8)
        ax_psd.set_xlabel("Hz", fontsize=8)
        ax_psd.set_yticks([])
        ax_psd.tick_params(labelsize=7)

    for row, area in enumerate(ROI):
        cc = c[c["area"] == area]
        for k in range(N_SHOW):
            ax_map, ax_psd = axes[row, 2 * k], axes[row, 2 * k + 1]
            if k >= len(cc):
                ax_map.axis("off"); ax_psd.axis("off")
                continue
            r = cc.iloc[k]
            comp = int(r["component"])
            title = f"{area} #{k + 1}: ICA{comp:03d}"
            if ref is not None:
                map_r = abs(np.corrcoef(ica.get_components()[:, comp], ref_map)[0, 1])
                env_r = np.corrcoef(band_envelope(ica, raw, comp, BAND[area]), ref_env)[0, 1]
                cand.loc[r.name, ["map_r", "env_r"]] = [map_r, env_r]
                title += f"\nvs current: map |r| {map_r:.2f}, envelope r {env_r:.2f}"
            draw(ax_map, ax_psd, ica, comp, r, area, title)
        if ref is not None:
            draw(axes[row, -2], axes[row, -1], ref_ica, ref[1], ref[2], area,
                 f"current: {ref[0]} ICA{ref[1]:03d}")
    what = "session ICA: top visual/alpha, bottom motor/beta" if ICA_NAME == "session" else \
        f"{ICA_NAME} localizer ICA" + (", last column = current session-ICA choice" if ref is not None else "")
    fig.suptitle(f"{sub}: candidates ({what}); dashed = 1/f line", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / f"candidates_{sub}{SUFFIX}.png", dpi=110)
    plt.close(fig)
    print("saved", OUT / f"candidates_{sub}{SUFFIX}.png", flush=True)

cand.round(3).to_csv(OUT / f"candidates_{SESSION}{SUFFIX}.csv", index=False)
print(cand.round(2).to_string(index=False))
