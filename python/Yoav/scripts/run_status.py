"""Show the progress of the latest execute_notebooks.py batch of one notebook.

Reads the status files <LOG_DIR>/<notebook>/status/<sub-NNN>.json and prints, per participant:
- waiting:  not started yet (the batch runs only --parallel N at a time);
- running:  which code cell (x of y) and section, and for how long;
- ok / failed: how long it took, and the error if it failed;
- DIED:     the status says running, but its process is gone (killed, or the job ended).
The last line says whether the batch is still running.

Usage (on UCloud, in the riseup_preproccesing folder; needs no uv environment):
    python3 python/Yoav/scripts/run_status.py 01_ica
"""

import json
import os
import sys
import time
from datetime import datetime
from pathlib import Path

LOG_DIR = Path(os.environ.get("RT_OUT_ROOT", "/work/YoavLevy#3923/derivatives/preproc")).parent.parent / "logs"


def alive(pid):
    """True if a process with this id exists on this machine."""
    if not pid:
        return False
    try:
        os.kill(pid, 0)  # signal 0: only checks that the process exists
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


def minutes_since(stamp):
    return (time.time() - datetime.strptime(stamp, "%Y-%m-%d %H:%M:%S").timestamp()) / 60


def main(notebook):
    folder = LOG_DIR / notebook / "status"
    batch_file = folder / "_batch.json"
    if not batch_file.exists():
        sys.exit(f"No batch of {notebook} found in {folder}")
    batch = json.loads(batch_file.read_text())
    print(f"{notebook}: batch started {batch['started']}, {len(batch['subjects'])} participants, "
          f"{batch['parallel']} at a time")

    counts = {}
    for subject in batch["subjects"]:
        path = folder / f"{subject}.json"
        s = json.loads(path.read_text()) if path.exists() else {"state": "waiting"}
        state = s["state"]
        if state == "running" and not alive(s.get("pid")):
            state = "DIED"
        counts[state] = counts.get(state, 0) + 1

        if state == "running":
            detail = (f"cell {s['cell']}/{s['n_cells']} ({s['section']}), "
                      f"{minutes_since(s['started']):.0f} min so far, this cell since {s['updated'][11:]}")
        elif state in ("ok", "failed"):
            detail = f"{s.get('minutes', '?')} min, finished {s['finished'][11:]}" + (
                f"\n{'':12}{s['error']}" if s.get("error") else "")
        elif state == "DIED":
            detail = f"stopped in cell {s['cell']}/{s['n_cells']} ({s['section']}), last update {s['updated']}"
        else:
            detail = ""
        print(f"  {subject:9} {state:8} {detail}")

    runner = "running" if alive(batch.get("runner_pid")) else "finished (or stopped)"
    summary = ", ".join(f"{n} {state}" for state, n in counts.items())
    print(f"Batch {runner}: {summary}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
