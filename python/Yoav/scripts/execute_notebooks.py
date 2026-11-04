"""Run a per-participant notebook (01_ica, 03_clean or 04_spatial_filter) for every participant
of a session, without opening it. Each run is saved, with its outputs, to
<LOG_DIR>/<notebook>/<notebook>_<sub-NNN>.ipynb; the notebook in the repo stays without outputs.

The participant is passed to the notebook through the environment variable RT_SUBJECT
(read in the notebook's "Participant" cell).

Usage (on UCloud, in the riseup_preproccesing folder):
    uv run python python/Yoav/scripts/execute_notebooks.py 01_ica session01
    uv run python python/Yoav/scripts/execute_notebooks.py 03_clean session01 sub-003   # one participant
    uv run python python/Yoav/scripts/execute_notebooks.py 02_corrmap session01         # all participants at once
"""

import os
import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parents[1]  # python/Yoav
sys.path.insert(0, str(HERE))
import config as cfg  # noqa: E402


def force_static_plots(nb):
    """Interactive plots (%matplotlib widget) do not show in a saved notebook: use static ones."""
    for cell in nb.cells:
        if cell.cell_type == "code":
            cell.source = cell.source.replace("%matplotlib widget", "%matplotlib inline")


def run(notebook, subject):
    nb = nbformat.read(HERE / f"{notebook}.ipynb", as_version=4)
    os.environ["RT_SUBJECT"] = subject  # the kernel inherits it
    out = cfg.LOG_DIR / notebook / f"{notebook}_{subject}.ipynb"
    out.parent.mkdir(parents=True, exist_ok=True)
    start = time.time()
    try:
        force_static_plots(nb)
        NotebookClient(nb, timeout=None, kernel_name="riseup", resources={"metadata": {"path": str(HERE)}}).execute()
        status = "ok"
    except Exception as error:  # keep going with the next participant, but save what ran
        status = f"FAILED: {type(error).__name__}: {str(error).splitlines()[-1] if str(error) else ''}"
    nbformat.write(nb, out)
    print(f"{notebook} {subject}: {status} ({(time.time() - start) / 60:.1f} min) -> {out}", flush=True)
    return status == "ok"


def main(notebook, group, subjects=None):
    if notebook == "02_corrmap":  # one run for all participants (the log is named after the session)
        subjects = [group]
    subjects = subjects or cfg.session_subjects(group)
    results = [run(notebook, s) for s in subjects]
    print(f"\n{sum(results)} of {len(results)} succeeded")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
