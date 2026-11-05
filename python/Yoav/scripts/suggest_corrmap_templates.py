"""Suggest CORRMAP templates (02_corrmap) with ICLabel, instead of choosing them by eye.

ICLabel (Pion-Tonachini et al., 2019; package mne-icalabel) gives every ICA component the
probability of being brain, muscle, eye, heart, line noise, channel noise or other. For each
artifact type, the suggested template is the component (over all participants with an ICA) with
the highest probability of that type:
- eye_blink:   eye, with Fp1 and Fp2 of the same sign and F7 and F8 of the same sign (symmetric map);
- lateral_eye: eye, with Fp1 and Fp2 and F7 and F8 of opposite signs (left/right map);
- heartbeat:   heart;
- muscle:      muscle.
A type is left out if no component reaches MIN_PROBABILITY: a weak template would make CORRMAP
remove components that are not that artifact.

Caveat: ICLabel was trained on average-referenced data; ours stay in the mastoid reference (as in
the reference pipeline), so mne-icalabel warns. The suggestions are a starting point: check the
template plots in 02_corrmap (step 2) and change CORRMAP_TEMPLATES in config.py if needed.

All probabilities are saved to <OUT_ROOT>/iclabel_probabilities.csv.

Usage (on UCloud, in the riseup_preproccesing folder):
    uv run python python/Yoav/scripts/suggest_corrmap_templates.py
"""

import sys
import warnings
from pathlib import Path

import mne
import pandas as pd
from mne_icalabel.iclabel import iclabel_label_components

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as cfg  # noqa: E402

MIN_PROBABILITY = 0.5
CLASSES = ["brain", "muscle", "eye", "heart", "line_noise", "channel_noise", "other"]  # ICLabel's order

mne.set_log_level("error")
warnings.filterwarnings("ignore", module="mne_icalabel")

rows = []
for ica_file in sorted(cfg.OUT_ROOT.glob("sub-*/sub-*_ica.fif")):
    subject = ica_file.name.removesuffix("_ica.fif")
    ica = mne.preprocessing.read_ica(ica_file)
    raw = mne.io.read_raw(ica_file.parent / f"{subject}_prep_raw.fif", preload=True)
    raw.filter(1, 100)  # ICLabel expects data high-passed at 1 Hz, like the data the ICA learned from
    probabilities = iclabel_label_components(raw, ica, inplace=False)  # (components, classes)
    maps = pd.DataFrame(ica.get_components(), index=ica.ch_names)  # (channels, components)
    for component, p in enumerate(probabilities):
        rows.append({"subject": subject, "component": component, **dict(zip(CLASSES, p.round(3))),
                     "same_sign_Fp": bool(maps.loc["Fp1", component] * maps.loc["Fp2", component] > 0),
                     "opposite_sign_F7_F8": bool(maps.loc["F7", component] * maps.loc["F8", component] < 0)})
    print(f"{subject}: {len(probabilities)} components labelled", flush=True)

table = pd.DataFrame(rows)
table.to_csv(cfg.OUT_ROOT / "iclabel_probabilities.csv", index=False)

candidates = {
    "eye_blink": ("eye", table["same_sign_Fp"] & ~table["opposite_sign_F7_F8"]),
    "lateral_eye": ("eye", ~table["same_sign_Fp"] & table["opposite_sign_F7_F8"]),
    "heartbeat": ("heart", True),
    "muscle": ("muscle", True),
}
print("\nSuggested CORRMAP_TEMPLATES (best component per type, probability in brackets):")
for label, (cls, condition) in candidates.items():
    pool = table[condition] if condition is not True else table
    best = pool.sort_values(cls, ascending=False).head(3)
    top = best.iloc[0]
    others = ", ".join(f"{r.subject} {r.component} ({r[cls]:.2f})" for _, r in best.iloc[1:].iterrows())
    if top[cls] >= MIN_PROBABILITY:
        print(f'    "{label}": ("{top.subject}", {top.component}),  # {cls} p = {top[cls]:.2f}; next: {others}')
    else:
        print(f"    # {label}: none (best {top.subject} {top.component}, {cls} p = {top[cls]:.2f} < {MIN_PROBABILITY})")
