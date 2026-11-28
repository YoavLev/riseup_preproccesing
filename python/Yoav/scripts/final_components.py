"""Figure of the chosen occipital alpha components (SPATIAL_FILTER_CHOICE) of one session, for presentations.

Two columns per participant (scalp map, spectrum), three participants per row. Each scalp map is seen from above
with the nose at the top, and labelled with its directions (front, back, L, R). Spectrum: the component on the
clean data (bad stretches left out), with the 1/f line (as component_scores.py) and the alpha band shaded; the
title gives the peak above the 1/f line. Participants without a choice are listed at the bottom (excluded).

Saves <OUT_ROOT>/spatial_filter_review/<session>/final_components_<session>.png

Usage (on UCloud, in riseup_preproccesing):
    uv run python python/Yoav/scripts/final_components.py session03
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mne
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as cfg  # noqa: E402
from component_scores import one_over_f_line  # noqa: E402

mne.set_log_level("error")
SESSION = sys.argv[1]
OUT = cfg.OUT_ROOT / "spatial_filter_review" / SESSION
OUT.mkdir(parents=True, exist_ok=True)
PER_ROW = 3
everyone = cfg.session_subjects(SESSION)
subjects = [s for s in everyone if s in cfg.SPATIAL_FILTER_CHOICE]
excluded = [s for s in everyone if s not in cfg.SPATIAL_FILTER_CHOICE]


def directions(ax):
    """Direction labels around a scalp map (MNE draws the head with the nose up)."""
    kw = dict(transform=ax.transAxes, ha="center", va="center", fontsize=8, color="0.25")
    ax.text(0.5, 1.07, "front", **kw)
    ax.text(0.5, -0.06, "back", **kw)
    ax.text(-0.06, 0.5, "L", fontweight="bold", **kw)
    ax.text(1.06, 0.5, "R", fontweight="bold", **kw)


n_rows = max(1, -(-len(subjects) // PER_ROW))
fig, axes = plt.subplots(n_rows, 2 * PER_ROW, figsize=(4.2 * PER_ROW, 2.8 * n_rows + 0.6),
                         gridspec_kw={"width_ratios": [1, 1.5] * PER_ROW}, squeeze=False)
for ax in axes.flat:
    ax.axis("off")
for k, sub in enumerate(subjects):
    ax_map, ax_psd = axes[k // PER_ROW, 2 * (k % PER_ROW)], axes[k // PER_ROW, 2 * (k % PER_ROW) + 1]
    ax_map.axis("on"); ax_psd.axis("on")
    comp = cfg.SPATIAL_FILTER_CHOICE[sub]
    raw = mne.io.read_raw(cfg.out_dir(sub) / f"{sub}_clean_raw.fif", preload=True)
    ica = mne.preprocessing.read_ica(cfg.out_dir(sub) / f"{sub}_visual_localizer_ica.fif")
    mne.viz.plot_topomap(ica.get_components()[:, comp], ica.info, axes=ax_map, show=False, contours=0, sensors=True)
    directions(ax_map)
    ax_map.set_title(f"{sub}\nvisual localizer ICA, component {comp}", fontsize=9, pad=14)
    psd = ica.get_sources(raw).pick([f"ICA{comp:03d}"]).compute_psd(fmin=1, fmax=40, picks="all", n_fft=1000,
                                                                    reject_by_annotation=True)
    f, p = psd.freqs, 10 * np.log10(psd.get_data()[0])
    line = one_over_f_line(f, p)
    lo, hi = cfg.ALPHA_BAND
    inside = (f >= lo) & (f <= hi)
    i = np.argmax((p - line)[inside])
    ax_psd.plot(f, p, color="k", lw=1.2)
    ax_psd.plot(f, line, color="0.6", lw=1, ls="--")
    ax_psd.axvspan(lo, hi, color="C0", alpha=0.15)
    ax_psd.set_title(f"peak {(p - line)[inside][i]:.1f} dB above 1/f at {f[inside][i]:g} Hz", fontsize=8)
    ax_psd.set_xlabel("frequency (Hz)", fontsize=8)
    ax_psd.set_yticks([])
    ax_psd.tick_params(labelsize=7)
    for s in ("top", "right"):
        ax_psd.spines[s].set_visible(False)
note = "Shaded: alpha band of the envelope; dashed: 1/f background. Maps seen from above, nose at the top."
if excluded:
    note += "  Excluded (no clear occipital alpha component): " + ", ".join(excluded) + "."
fig.suptitle(f"{SESSION}: chosen occipital alpha components", fontsize=12)
fig.text(0.5, 0.005, note, ha="center", fontsize=8.5, color="0.3")
fig.tight_layout(rect=(0, 0.02, 1, 0.98))
fig.savefig(OUT / f"final_components_{SESSION}.png", dpi=120)
print("saved", OUT / f"final_components_{SESSION}.png")
