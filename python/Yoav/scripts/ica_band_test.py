"""Test: does the ICA give cleaner components when it learns from 1-100 Hz data at 500 Hz?

01_ica trains the ICA on a copy filtered to 1-40 Hz and downsampled to 250 Hz (the reference
pipeline, Stropahl et al. 2018). The lab's s02_erp_ica notebook (RiseTogether-analysis) trains
on 1-100 Hz at 500 Hz with a 50/100 Hz notch, and its components look cleaner. This script refits
one participant with that band and changes nothing else:
- same input: <sub>_prep_raw.fif from 01_ica (bad channels interpolated, average reference,
  breaks and BAD_amplitude stretches marked);
- same filters (EEGLAB-style FIR, Hann window). The filter orders are doubled because the
  sampling rate is doubled, so the transition bands are as wide in Hz as in 01_ica;
- same 1 s dummy epochs, joint-probability rejection at 3 SD, extended Infomax, number of
  components and random seed.
Only new: the 50 and 100 Hz notch (mains noise is inside the band now), the 100 Hz low-pass and
no downsampling.

To compare the two fits, every component of both gets an ICLabel label (computed on the band that
ICA learned from), and its best match in the other fit (largest |correlation| of the scalp maps).

Outputs, in derivatives/ica_band_test/<sub>/ (nothing in derivatives/preproc is changed):
- <sub>_ica_1-100Hz.fif                the new ICA;
- <sub>_maps_1-40Hz.png, _1-100Hz.png  all maps of each fit: ICLabel class and probability, and the
                                       best match in the other fit; red = removed by 02_corrmap
                                       (1-40 Hz fit only);
- <sub>_band_comparison.csv            one row per component of each fit (labels, matches), plus
                                       the fit numbers (rejected epochs, iterations).

Usage (on UCloud, in the riseup_preproccesing folder):
    uv run python python/Yoav/scripts/ica_band_test.py sub-016
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

LP_WIDE = 100  # low-pass of the new training copy [Hz]
NOTCH = [50, 100]  # mains frequency and its harmonic [Hz]
CLASSES = ["brain", "muscle", "eye", "heart", "line_noise", "channel_noise", "other"]  # ICLabel's order

mne.set_log_level("error")
warnings.filterwarnings("ignore", module="mne_icalabel")

subject = sys.argv[1]  # e.g. sub-016
src_dir = cfg.out_dir(subject)
out = cfg.OUT_ROOT.parent / "ica_band_test" / subject
out.mkdir(parents=True, exist_ok=True)

raw = mne.io.read_raw(src_dir / f"{subject}_prep_raw.fif", preload=True)
ica_old = mne.preprocessing.read_ica(src_dir / f"{subject}_ica.fif")
removed_file = src_dir / f"{subject}_badcomps.csv"
removed = pd.read_csv(removed_file)["component"].astype(int).tolist() if removed_file.exists() else []
print(f"{subject}: {raw.info['sfreq']} Hz, {raw.times[-1]:.0f} s; 1-40 Hz fit: {ica_old.n_components_} components",
      flush=True)

# --- 1. The new training copy: notch, then 1-100 Hz, at 500 Hz ---
raw_wide = raw.copy()
raw_wide.notch_filter(NOTCH, picks=["eeg", "eog"])
scale = raw.info["sfreq"] / cfg.ICA_SRATE  # 2: twice the samples per second, so twice the filter order
firws_filter(raw_wide, LP_WIDE, "lowpass", int(cfg.ICA_LP_ORDER * scale))
firws_filter(raw_wide, cfg.ICA_HP, "highpass", int(cfg.ICA_HP_ORDER * scale))

# --- 2. Dummy epochs and rejection, as in 01_ica ---
epochs = mne.make_fixed_length_epochs(raw_wide, duration=cfg.ICA_EPOCH_LENGTH, reject_by_annotation=True,
                                      preload=True)
keep = jointprob_keep(epochs.get_data(picks="eeg"), cfg.ICA_PRUNE)
epochs_ica = epochs[np.where(keep)[0]]
print(f"Dummy epochs: {len(epochs)}, rejected {np.sum(~keep)}", flush=True)

# --- 3. ICA, as in 01_ica (same number of components and seed) ---
ica_new = ICA(n_components=ica_old.n_components_, method="infomax", fit_params=dict(extended=True),
              max_iter="auto", random_state=cfg.RANDOM_STATE)
ica_new.fit(epochs_ica, picks="eeg")
ica_new.save(out / f"{subject}_ica_1-100Hz.fif", overwrite=True)
print(f"1-100 Hz fit: {ica_new.n_components_} components, {ica_new.n_iter_} iterations", flush=True)


# --- 4. ICLabel, on the band each ICA learned from (1 s pieces without the BAD_ stretches) ---
def iclabel(ica, data):
    pieces = mne.make_fixed_length_epochs(data, duration=1.0, reject_by_annotation=True, preload=True)
    p = iclabel_label_components(pieces, ica, inplace=False)  # (components, classes)
    return pd.DataFrame(p, columns=CLASSES)


raw_narrow = raw.copy().filter(1, cfg.ICA_LP)  # as in suggest_corrmap_templates.py
labels = {"1-40Hz": iclabel(ica_old, raw_narrow), "1-100Hz": iclabel(ica_new, raw_wide)}

# --- 5. Best match of every component in the other fit (|correlation| of the scalp maps) ---
maps = {"1-40Hz": ica_old.get_components(), "1-100Hz": ica_new.get_components()}  # (channels, components)
n_old = maps["1-40Hz"].shape[1]
r = np.abs(np.corrcoef(maps["1-40Hz"].T, maps["1-100Hz"].T)[:n_old, n_old:])  # (old, new)

rows = []
for fit, other, rr in [("1-40Hz", "1-100Hz", r), ("1-100Hz", "1-40Hz", r.T)]:
    lab = labels[fit]
    for c in range(len(lab)):
        best = int(rr[c].argmax())
        rows.append({"fit": fit, "component": c, "label": lab.iloc[c].idxmax(), "p": round(lab.iloc[c].max(), 2),
                     **lab.iloc[c].round(3).to_dict(), "match_in_" + other.replace("-", "_"): best,
                     "match_r": round(rr[c, best], 2), "removed_by_corrmap": fit == "1-40Hz" and c in removed})
table = pd.DataFrame(rows)
fit_info = pd.DataFrame([
    {"fit": "1-40Hz", "component": "fit", "n_iter": ica_old.n_iter_},
    {"fit": "1-100Hz", "component": "fit", "n_iter": ica_new.n_iter_, "n_dummy_epochs": len(epochs),
     "n_rejected_epochs": int(np.sum(~keep))},
])
pd.concat([table, fit_info]).to_csv(out / f"{subject}_band_comparison.csv", index=False)

# --- 6. One PNG per fit with all maps ---
for fit, other, ica in [("1-40Hz", "1_100Hz", ica_old), ("1-100Hz", "1_40Hz", ica_new)]:
    t = table[table["fit"] == fit]
    notes = {row.component: f"{row.label} {row.p:.2f}\n≈{int(getattr(row, 'match_in_' + other)):03d} ({row.match_r:.2f})"
             for row in t.itertuples()}
    title = (f"{subject}, ICA trained on {fit}: ICLabel class and p; ≈ best match in the other fit (|r|)"
             + ("; red = removed by 02_corrmap" if fit == "1-40Hz" else ""))
    figure = plot_component_maps(ica, title, notes=notes, removed=removed if fit == "1-40Hz" else ())
    figure.savefig(out / f"{subject}_maps_{fit}.png", dpi=100)
    plt.close(figure)

# --- 7. Summary for the log ---
print("\nICLabel classes (most likely class per component):")
print(table.pivot_table(index="label", columns="fit", values="component", aggfunc="count", fill_value=0))
print("\nComponents with a close match in the other fit (|r| >= 0.9):")
print(table.groupby("fit")["match_r"].apply(lambda m: int((m >= 0.9).sum())))
print("Saved to", out)
