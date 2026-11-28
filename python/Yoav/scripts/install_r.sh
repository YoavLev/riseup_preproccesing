#!/usr/bin/env bash
# Installs R and the R packages for 06_crqa.R and 07_crqa_plots.R on a UCloud job. Run once per job (after setup_job.sh):
#   bash "/work/YoavLevy#3923/riseup_preproccesing/python/Yoav/scripts/install_r.sh"
#
# R itself (apt, Ubuntu's r-base 4.3.3) is part of the job's system and is gone after every job, so it is
# installed again each time (~1 min). The R packages (crqa, data.table and their dependencies) go to a
# library in /work, which is kept: they are compiled only the first time (~5 min), later runs only check them.
set -euo pipefail

WORK="/work/YoavLevy#3923"
RLIB="$WORK/R-lib"          # R package library, kept between jobs
mkdir -p "$RLIB" "$WORK/logs/run"

# 1. R and the system libraries the packages need to compile (sudo works on UCloud jobs)
if ! command -v Rscript >/dev/null; then
  sudo DEBIAN_FRONTEND=noninteractive apt-get update -qq
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
    r-base-core r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev libglpk-dev libgmp-dev > /dev/null
fi

# 2. Tell R where the package library is, in every new shell of this job (the home is wiped with the job)
echo "R_LIBS_USER=$RLIB" > "$HOME/.Renviron"
export R_LIBS_USER="$RLIB"

# 3. The R packages, only those that are missing
Rscript -e '
options(repos = "https://cloud.r-project.org", Ncpus = 8)
lib <- Sys.getenv("R_LIBS_USER")
for (p in c("data.table", "crqa", "lme4")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, lib = lib)
}
library(crqa)
cat(R.version.string, "| crqa", format(packageVersion("crqa")), "| data.table", format(packageVersion("data.table")), "\n")
'
echo "R ready. Run: Rscript 06_crqa.R session02,session03,... [cores], then Rscript 07_crqa_plots.R"
