"""Test: does a change to the ICA give cleaner components? Refits one participant with one change.

The lab's s02_erp_ica notebook (RiseTogether-analysis) gets cleaner-looking components than
01_ica. This script refits one participant with one of its settings and changes nothing else, so
the effect of that one setting can be seen. Variants:
- "1-100Hz": train on 1-100 Hz at 500 Hz with a 50/100 Hz notch (s02), instead of 1-40 Hz at
  250 Hz (01_ica, from the reference pipeline). The filter orders are doubled because the sampling
  rate is doubled, so the transition bands are as wide in Hz as in 01_ica. Still extended Infomax.
- "picard": the same training data as 01_ica (1-40 Hz, 250 Hz), but fitted with Picard
  (extended, ortho=False, max_iter=1000, as in s02) instead of MNE's Infomax. Picard fits the
  same model (extended Infomax) but stops when it has converged. MNE's Infomax stops at 500
  iterations, which every 01_ica fit reached, so they may not have converged.
Everything else is as in 01_ica: input <sub>_prep_raw.fif (bad channels interpolated, average
reference, breaks and BAD_amplitude marked), EEGLAB-style FIR filters, 1 s dummy epochs,
joint-probability rejection at 3 SD, same number of components, same random seed.

To compare with the current fit (01_ica), every component of both fits gets an ICLabel label and
its best match in the other fit (largest |correlation| of the scalp maps). ICLabel gets the band
the ICA learned from (1-40 Hz for the current and the picard fit, 1-100 Hz for the 1-100Hz fit).

Outputs, in derivatives/ica_test/<sub>/<variant>/ (nothing in derivatives/preproc is changed):
- <sub>_ica_<variant>.fif      the new ICA;
- <sub>_maps_current.png       all maps of the current fit: ICLabel class and p, best match in the
                               new fit (|r|); red = removed by 02_corrmap;
- <sub>_maps_<variant>.png     the same for the new fit (match in the current fit);
- <sub>_comparison.csv         one row per component of each fit (labels, matches), plus the fit
                               numbers (rejected epochs, iterations).

Usage (on UCloud, in the riseup_preproccesing folder):
    uv run python python/Yoav/scripts/ica_test.py sub-016 picard
"""

import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # no screen: only save the figures
import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd
from mne.preprocessing import ICA
from mne_icalabel.iclabel import iclabel_label_components

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as cfg  # noqa: E402
from utils import firws_filter, jointprob_keep, plot_component_maps  # noqa: E402

VARIANTS = ["1-100Hz", "picard"]
LP_WIDE = 100  # low-pass of the 1-100Hz training copy [Hz]
NOTCH = [50, 100]  # mains frequency and its harmonic [Hz]
PICARD_MAX_ITER = 1000  # as in s02_erp_ica
CLASSES = ["brain", "muscle", "eye", "heart", "line_noise", "channel_noise", "other"]  # ICLabel's order

mne.set_log_level("error")
warnings.filterwarnings("ignore", module="mne_icalabel")

subject, variant = sys.argv[1], sys.argv[2]  # e.g. sub-016 picard
assert variant in VARIANTS, f"variant must be one of {VARIANTS}"
src_dir = cfg.out_dir(subject)
out = cfg.OUT_ROOT.parent / "ica_test" / subject / variant
out.mkdir(parents=True, exist_ok=True)

raw = mne.io.read_raw(src_dir / f"{subject}_prep_raw.fif", preload=True)
ica_old = mne.preprocessing.read_ica(src_dir / f"{subject}_ica.fif")
removed_file = src_dir / f"{subject}_badcomps.csv"
removed = pd.read_csv(removed_file)["component"].astype(int).tolist() if removed_file.exists() else []
print(f"{subject}, variant {variant}: {raw.info['sfreq']} Hz, {raw.times[-1]:.0f} s; "
      f"current fit: {ica_old.n_components_} components", flush=True)

# --- 1. The training copy ---
train = raw.copy()
if variant == "1-100Hz":  # notch, then 1-100 Hz, at 500 Hz
    train.notch_filter(NOTCH, picks=["eeg", "eog"])
    scale = raw.info["sfreq"] / cfg.ICA_SRATE  # 2: twice the samples per second, so twice the filter order
    firws_filter(train, LP_WIDE, "lowpass", int(cfg.ICA_LP_ORDER * scale))
    firws_filter(train, cfg.ICA_HP, "highpass", int(cfg.ICA_HP_ORDER * scale))
else:  # as in 01_ica: low-pass 40 Hz, 250 Hz, high-pass 1 Hz
    firws_filter(train, cfg.ICA_LP, "lowpass", cfg.ICA_LP_ORDER)
    train.resample(cfg.ICA_SRATE)
    firws_filter(train, cfg.ICA_HP, "highpass", cfg.ICA_HP_ORDER)

# --- 2. Dummy epochs and rejection, as in 01_ica ---
epochs = mne.make_fixed_length_epochs(train, duration=cfg.ICA_EPOCH_LENGTH, reject_by_annotation=True,
                                      preload=True)
keep = jointprob_keep(epochs.get_data(picks="eeg"), cfg.ICA_PRUNE)
epochs_ica = epochs[np.where(keep)[0]]
print(f"Dummy epochs: {len(epochs)}, rejected {np.sum(~keep)}", flush=True)

# --- 3. ICA (same number of components and seed as 01_ica) ---
if variant == "picard":
    ica_new = ICA(n_components=ica_old.n_components_, method="picard", max_iter=PICARD_MAX_ITER,
                  fit_params=dict(extended=True, ortho=False), random_state=cfg.RANDOM_STATE)
else:
    ica_new = ICA(n_components=ica_old.n_components_, method="infomax", fit_params=dict(extended=True),
                  max_iter="auto", random_state=cfg.RANDOM_STATE)
ica_new.fit(epochs_ica, picks="eeg")
ica_new.save(out / f"{subject}_ica_{variant}.fif", overwrite=True)
limit = PICARD_MAX_ITER if variant == "picard" else 500  # MNE's "auto" for Infomax
print(f"{variant} fit: {ica_new.n_components_} components, {ica_new.n_iter_} iterations "
      f"({'stopped at the limit' if ica_new.n_iter_ >= limit else 'converged before the limit'})", flush=True)


# --- 4. ICLabel, on the band each ICA learned from (1 s pieces without the BAD_ stretches) ---
def iclabel(ica, data):
    pieces = mne.make_fixed_length_epochs(data, duration=1.0, reject_by_annotation=True, preload=True)
    p = iclabel_label_components(pieces, ica, inplace=False)  # (components, classes)
    return pd.DataFrame(p, columns=CLASSES)


raw_narrow = raw.copy().filter(1, cfg.ICA_LP)  # as in suggest_corrmap_templates.py
labels = {"current": iclabel(ica_old, raw_narrow),
          variant: iclabel(ica_new, train if variant == "1-100Hz" else raw_narrow)}

# --- 5. Best match of every component in the other fit (|correlation| of the scalp maps) ---
maps_old, maps_new = ica_old.get_components(), ica_new.get_components()  # (channels, components)
n_old = maps_old.shape[1]
r = np.abs(np.corrcoef(maps_old.T, maps_new.T)[:n_old, n_old:])  # (current, new)

rows = []
for fit, rr in [("current", r), (variant, r.T)]:
    lab = labels[fit]
    for c in range(len(lab)):
        best = int(rr[c].argmax())
        rows.append({"fit": fit, "component": c, "label": lab.iloc[c].idxmax(), "p": round(lab.iloc[c].max(), 2),
                     **lab.iloc[c].round(3).to_dict(), "match_in_other_fit": best, "match_r": round(rr[c, best], 2),
                     "removed_by_corrmap": fit == "current" and c in removed})
table = pd.DataFrame(rows)
fit_info = pd.DataFrame([
    {"fit": "current", "component": "fit", "n_iter": ica_old.n_iter_},
    {"fit": variant, "component": "fit", "n_iter": ica_new.n_iter_, "n_dummy_epochs": len(epochs),
     "n_rejected_epochs": int(np.sum(~keep))},
])
pd.concat([table, fit_info]).to_csv(out / f"{subject}_comparison.csv", index=False)

# --- 6. One PNG per fit with all maps ---
for fit, ica in [("current", ica_old), (variant, ica_new)]:
    t = table[table["fit"] == fit]
    notes = {row.component: f"{row.label} {row.p:.2f}\n≈{row.match_in_other_fit:03d} ({row.match_r:.2f})"
             for row in t.itertuples()}
    other = variant if fit == "current" else "current"
    title = (f"{subject}, {fit} fit: ICLabel class and p; ≈ best match in the {other} fit (|r|)"
             + ("; red = removed by 02_corrmap" if fit == "current" else ""))
    figure = plot_component_maps(ica, title, notes=notes, removed=removed if fit == "current" else ())
    figure.savefig(out / f"{subject}_maps_{fit}.png", dpi=100)
    plt.close(figure)

# --- 7. Summary for the log ---
print("\nICLabel classes (most likely class per component):")
print(table.pivot_table(index="label", columns="fit", values="component", aggfunc="count", fill_value=0))
print("\nComponents with a close match in the other fit (|r| >= 0.9):")
print(table.groupby("fit")["match_r"].apply(lambda m: int((m >= 0.9).sum())))
print("Saved to", out)
