#!/usr/bin/env bash
# Run once at the start of every UCloud job (the job's home is wiped, /work is kept):
#   bash "/work/YoavLevy#3923/riseup_preproccesing/python/Yoav/scripts/setup_job.sh"
# Installs uv and the Python environment of riseup_preproccesing, and registers it as a
# Jupyter kernel for VS Code. Mount the dataset (/work/RiseTogetherEEG or /work/VideoGameParadigm-EEG) and /work/YoavLevy#3923 in the job.
set -euo pipefail

WORK="/work/YoavLevy#3923"
REPO="$WORK/riseup_preproccesing"
# uv cannot handle '#' in some paths: keep the environment outside the '#' folder (in the job's home;
# it is rebuilt each job, which takes about a minute from uv's cache)
export UV_PROJECT_ENVIRONMENT="$HOME/riseup-venv"
export UV_CACHE_DIR="$WORK/.uv-cache"

command -v uv >/dev/null || curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="$HOME/.local/bin:$PATH"

cd "$REPO"
uv sync
uv run python -m ipykernel install --user --name riseup --display-name "riseup (uv)"
# interactive plots in JupyterLab (%matplotlib widget): ipympl must also be in JupyterLab's own Python
if [ -x /opt/conda/bin/pip ]; then /opt/conda/bin/pip install -q ipympl; fi

mkdir -p "$WORK/derivatives/preproc" "$WORK/logs"
test -d /work/RiseTogetherEEG/VideoGameParadigm-EEG || test -d /work/VideoGameParadigm-EEG || echo "WARNING: the BIDS dataset is not mounted"
echo "Ready. In new shells run:  export UV_PROJECT_ENVIRONMENT=\"$UV_PROJECT_ENVIRONMENT\" PATH=\"\$HOME/.local/bin:\$PATH\""
