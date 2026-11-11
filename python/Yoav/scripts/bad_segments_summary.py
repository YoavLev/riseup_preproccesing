"""How much each participant loses to BAD_amplitude (01_ica step 3b), from the saved <sub>_bad_segments.csv.

Usage (from riseup_preproccesing/):  uv run python python/Yoav/scripts/bad_segments_summary.py session03
Prints one row per participant (seconds and % of trial time lost, worst trial) and the trials losing more than 10%.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # python/Yoav, for config
import config as cfg

session = sys.argv[1]
tables = []
for subject in cfg.session_subjects(session):
    path = cfg.out_dir(subject) / f"{subject}_bad_segments.csv"
    if path.exists():
        tables.append(pd.read_csv(path))
    else:
        print(f"{subject}: no {path.name} yet (run 01_ica first)")
if not tables:
    sys.exit()
trials = pd.concat(tables)

per_subject = trials.groupby("subject").agg(trial_s=("trial_s", "sum"), bad_s=("bad_s", "sum"),
                                            worst_trial_percent=("bad_percent", "max"))
per_subject["bad_percent"] = (100 * per_subject["bad_s"] / per_subject["trial_s"]).round(1)
print(f"\n{session}: time marked BAD_amplitude, per participant")
print(per_subject[["trial_s", "bad_s", "bad_percent", "worst_trial_percent"]].to_string())

big = trials[trials["bad_percent"] > 10]
print(f"\nTrials losing more than 10%: {len(big)}")
if len(big):
    print(big.to_string(index=False))
