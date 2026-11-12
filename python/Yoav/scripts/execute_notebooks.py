"""Run a per-participant notebook (01_ica, 03_clean or 04_spatial_filter) for every participant
of a session, without opening it. Each run is saved, with its outputs, to
<LOG_DIR>/<notebook>/<notebook>_<sub-NNN>.ipynb; the notebook in the repo stays without outputs.

The participant is passed to the notebook through the environment variable RT_SUBJECT
(read in the notebook's "Participant" cell).

--parallel N runs N participants at the same time (default 1: one after the other).

Progress: every run keeps a small status file, <LOG_DIR>/<notebook>/status/<sub-NNN>.json
(state, process id, which cell is running, start and end time). It is updated when each cell
starts, so it is current even while a long cell runs. Read them with run_status.py:
    uv run python python/Yoav/scripts/run_status.py 01_ica

Usage (on UCloud, in the riseup_preproccesing folder; with nohup so the run survives a lost connection):
    nohup uv run python python/Yoav/scripts/execute_notebooks.py 01_ica session03 --parallel 3 \
        > "/work/YoavLevy#3923/logs/run/01_ica_session03.txt" 2>&1 &
    uv run python python/Yoav/scripts/execute_notebooks.py 03_clean session01 sub-003   # one participant
    uv run python python/Yoav/scripts/execute_notebooks.py 02_corrmap session01         # all participants at once
"""

import json
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import nbformat
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parents[1]  # python/Yoav
sys.path.insert(0, str(HERE))
import config as cfg  # noqa: E402


def status_dir(notebook):
    return cfg.LOG_DIR / notebook / "status"


def write_status(notebook, subject, **fields):
    """Update the run's status file (written to a temporary file first, so a reader never sees half a file)."""
    path = status_dir(notebook) / f"{subject}.json"
    status = json.loads(path.read_text()) if path.exists() else {}
    status.update(fields, updated=time.strftime("%Y-%m-%d %H:%M:%S"))
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(status, indent=1))
    tmp.replace(path)


def section_names(nb):
    """For each cell: the last markdown heading above it (e.g. '6. ICA training copy'), to show progress."""
    names, current = [], ""
    for cell in nb.cells:
        if cell.cell_type == "markdown" and cell.source.lstrip().startswith("#"):
            current = cell.source.lstrip().splitlines()[0].lstrip("#").strip()
        names.append(current)
    return names


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
    sections = section_names(nb)
    code_cells = [i for i, c in enumerate(nb.cells) if c.cell_type == "code"]
    write_status(notebook, subject, state="running", pid=os.getpid(), started=time.strftime("%Y-%m-%d %H:%M:%S"),
                 finished=None, cell=0, n_cells=len(code_cells), section="", error=None)

    def on_cell_start(cell, cell_index):
        # called for every cell (also markdown); progress = number of code cells up to this one.
        # Never let a progress update stop the notebook.
        try:
            done = sum(1 for i in code_cells if i <= cell_index)
            write_status(notebook, subject, cell=done, section=sections[cell_index])
        except Exception:
            pass

    start = time.time()
    error = None
    try:
        force_static_plots(nb)
        NotebookClient(nb, timeout=None, kernel_name="riseup", resources={"metadata": {"path": str(HERE)}},
                       on_cell_start=on_cell_start).execute()
    except Exception as e:  # keep going with the next participant, but save what ran
        last_line = re.sub(r"\x1b\[[0-9;]*m", "", str(e)).strip().splitlines()[-1] if str(e).strip() else ""  # no colour codes
        error = f"{type(e).__name__}: {last_line}"
    nbformat.write(nb, out)
    minutes = (time.time() - start) / 60
    write_status(notebook, subject, state="failed" if error else "ok", error=error, minutes=round(minutes, 1),
                 finished=time.strftime("%Y-%m-%d %H:%M:%S"))
    print(f"{notebook} {subject}: {'FAILED: ' + error if error else 'ok'} ({minutes:.1f} min) -> {out}", flush=True)
    return error is None


def main(notebook, group, subjects=None, parallel=1):
    if notebook == "02_corrmap":  # one run for all participants (the log is named after the session)
        subjects = [group]
    subjects = subjects or cfg.session_subjects(group)
    # the list of this batch, so run_status.py can show who is still waiting
    status_dir(notebook).mkdir(parents=True, exist_ok=True)
    for s in subjects:
        (status_dir(notebook) / f"{s}.json").unlink(missing_ok=True)  # status of an earlier run
        write_status(notebook, s, state="waiting", pid=None, cell=0, section="", error=None)
    (status_dir(notebook) / "_batch.json").write_text(json.dumps(
        {"subjects": subjects, "parallel": parallel, "runner_pid": os.getpid(),
         "started": time.strftime("%Y-%m-%d %H:%M:%S")}, indent=1))

    with ProcessPoolExecutor(max_workers=parallel) as pool:
        results = list(pool.map(run, [notebook] * len(subjects), subjects))
    print(f"\n{sum(results)} of {len(results)} succeeded", flush=True)


if __name__ == "__main__":
    args = sys.argv[1:]
    parallel = 1
    if "--parallel" in args:
        i = args.index("--parallel")
        parallel = int(args[i + 1])
        del args[i:i + 2]
    if len(args) < 2:
        sys.exit(__doc__)
    main(args[0], args[1], args[2:], parallel)
