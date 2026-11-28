"""Shortlist of occipital alpha components of the visual localizer ICA, for choosing SPATIAL_FILTER_CHOICE by eye.

Input: the numbers of scripts/component_scores.py <session>. Rule (constants in config.py): a component is a
candidate if its strongest electrode is occipital (OCCIPITAL: O1, O2), at least POSTERIOR_MIN_SHARE (35 %) of its
map's squared weights lie in the posterior electrodes (POSTERIOR), its spectrum is brain-like (20-40 Hz slope below
MAX_SLOPE_20_40, -0.5; muscle is flat) and its alpha peak is at least MIN_ALPHA_DB (1 dB) above the 1/f background.
Candidates are ranked by the alpha peak. The top N_SHOW are drawn (scalp map + spectrum with the 1/f line and the
alpha band) so the final choice is made by eye. A participant without a candidate is listed at the end (excluded
unless a component is chosen by eye anyway, with a reason).

Saves <OUT_ROOT>/spatial_filter_review/<session>/candidates_<sub>.png and candidates_<session>.csv.

Usage (on UCloud, in riseup_preproccesing, after component_scores.py):
    uv run python python/Yoav/scripts/component_candidates.py session03
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
from component_scores import one_over_f_line  # noqa: E402

mne.set_log_level("error")
SESSION = sys.argv[1]
OUT = cfg.OUT_ROOT / "spatial_filter_review" / SESSION
OUT.mkdir(parents=True, exist_ok=True)
N_SHOW = 3

scores = pd.read_csv(cfg.OUT_ROOT / "component_scores" / f"{SESSION}_visual_ica.csv")
scores["top1"] = scores["top3"].str.split().str[0]
ok = (scores["top1"].isin(cfg.OCCIPITAL) & (scores["posterior_share"] >= cfg.POSTERIOR_MIN_SHARE)
      & (scores["slope_20_40"] < cfg.MAX_SLOPE_20_40) & (scores["alpha_db"] >= cfg.MIN_ALPHA_DB))
cand = scores[ok].sort_values(["subject", "alpha_db"], ascending=[True, False])
cand = cand.assign(rank=cand.groupby("subject").cumcount() + 1)
cand = cand[cand["rank"] <= N_SHOW].drop(columns="top1")

for sub, c in cand.groupby("subject"):
    raw = mne.io.read_raw(cfg.out_dir(sub) / f"{sub}_clean_raw.fif", preload=True)
    ica = mne.preprocessing.read_ica(cfg.out_dir(sub) / f"{sub}_visual_localizer_ica.fif")
    fig, axes = plt.subplots(1, 2 * N_SHOW, figsize=(4.2 * N_SHOW, 3.4),
                             gridspec_kw={"width_ratios": [1, 1.4] * N_SHOW}, squeeze=False)
    for k in range(N_SHOW):
        ax_map, ax_psd = axes[0, 2 * k], axes[0, 2 * k + 1]
        if k >= len(c):
            ax_map.axis("off"); ax_psd.axis("off")
            continue
        r = c.iloc[k]
        comp = int(r["component"])
        mne.viz.plot_topomap(ica.get_components()[:, comp], ica.info, axes=ax_map, show=False, contours=0)
        ax_map.set_title(f"#{k + 1}: ICA{comp:03d}\n{r['top3']}", fontsize=8)
        psd = ica.get_sources(raw).pick([f"ICA{comp:03d}"]).compute_psd(fmin=1, fmax=40, picks="all", n_fft=1000,
                                                                        reject_by_annotation=True)
        f, p = psd.freqs, 10 * np.log10(psd.get_data()[0])
        ax_psd.plot(f, p, color="k", lw=1)
        ax_psd.plot(f, one_over_f_line(f, p), color="0.6", lw=1, ls="--")
        ax_psd.axvspan(*cfg.ALPHA_BAND, color="C0", alpha=0.15)
        ax_psd.set_title(f"alpha peak {r['alpha_db']:.1f} dB at {r['alpha_hz']:g} Hz; slope {r['slope_20_40']:.1f}",
                         fontsize=8)
        ax_psd.set_xlabel("Hz", fontsize=8)
        ax_psd.set_yticks([])
        ax_psd.tick_params(labelsize=7)
    fig.suptitle(f"{sub}: occipital alpha candidates (visual localizer ICA); dashed = 1/f line", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / f"candidates_{sub}.png", dpi=110)
    plt.close(fig)
    print("saved", OUT / f"candidates_{sub}.png", flush=True)

cand.round(3).to_csv(OUT / f"candidates_{SESSION}.csv", index=False)
print(cand.round(2).to_string(index=False))
none = sorted(set(cfg.session_subjects(SESSION)) - set(cand["subject"]))
print("no candidate (check by eye, otherwise excluded):", ", ".join(none) if none else "none")
