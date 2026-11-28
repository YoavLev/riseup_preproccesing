# 06_crqa.R: cRQA of all sessions with the settings decided on 2026-11-25 (after a parameter exploration).
#
# Input: the window folders of 05_windows with the decided smoothing (config.py: CRQA_SMOOTH_HZ = 1 Hz, zero-phase
# FIR, then CRQA_SFREQ = 5 Hz): derivatives/crqa/<session>_smooth1hz_5hz/ (windows.csv, envelopes.csv).
#
# Settings (one set for every session, pair, role and window):
#   - each envelope z-scored per window;
#   - embedding m = 1 (no embedding): the false-nearest-neighbour share of the envelopes stayed at 80-90 % for every
#     dimension (they behave like noise at this time scale), so no embedding dimension is supported by the data.
#     With m = 1 the delay plays no role;
#   - radius: one per area, the median over all real windows (all sessions, all roles) of the 5 % quantile of the
#     distances between A's and B's values (target recurrence ~5 %, the upper end of the range of Coco et al., 2021);
#     fixed, so recurrence can differ between roles;
#   - minimum diagonal and vertical line: 0.4 s (2 samples at 5 Hz);
#   - primary measure: lagRR_2s = %REC within +-2 s of the same moment (diagonals |lag| <= 2 s). Control: RR, the
#     recurrence rate of the whole plot (any lag), and DET.
#   - surrogate pairs (joint trials): A's envelope in one joint trial with B's envelope in a joint trial of another
#     block of the same session (same trial number), cut to the shorter. Same people, same task, not the same moment.
#
# Usage (on UCloud, after scripts/install_r.sh):  Rscript 06_crqa.R session02,session03,session04,session05,session06 [cores]
# Output: derivatives/crqa/final_smooth1hz_5hz/  parameters.csv, windows.csv, results.csv, profiles.csv

suppressMessages({
  library(crqa)
  library(data.table)
  library(parallel)
  library(Matrix)
})
setDTthreads(1)   # data.table's OpenMP threads make forked mclapply() workers hang

args <- commandArgs(trailingOnly = TRUE)
sessions <- strsplit(if (length(args) >= 1) args[1] else "session03,session04", ",")[[1]]
cores <- if (length(args) >= 2) as.integer(args[2]) else 8
suffix <- Sys.getenv("RT_CRQA_SUFFIX", "smooth1hz_5hz")
root <- Sys.getenv("RT_CRQA_ROOT", "/work/YoavLevy#3923/derivatives/crqa")
out <- file.path(root, paste0("final_", suffix))
dir.create(out, showWarnings = FALSE)
RNGkind("L'Ecuyer-CMRG")
set.seed(97)

# --- Settings -------------------------------------------------------------------------------------------------
embed <- 1
delay <- 1            # no role with embed = 1
target <- 5           # target recurrence rate [%] for the radius
line_s <- 0.4         # minimum diagonal / vertical line [s]
band_s <- 2           # lagRR: recurrence within +-2 s of the same moment
profile_s <- 10       # diagonal profile saved for -10..10 s

# --- Data -----------------------------------------------------------------------------------------------------
read_session <- function(s) {
  dir <- file.path(root, paste0(s, "_", suffix))
  w <- fread(file.path(dir, "windows.csv"))
  e <- fread(file.path(dir, "envelopes.csv"))
  w[, uid := paste(s, window_id, sep = ":")]   # window ids are only unique within a session
  e[, uid := paste(s, window_id, sep = ":")]
  list(w = w, e = e)
}
data <- lapply(sessions, read_session)
windows <- rbindlist(lapply(data, `[[`, "w"), fill = TRUE)
envelopes <- rbindlist(lapply(data, `[[`, "e"), fill = TRUE)
rm(data)
envelopes[, x_a := as.numeric(scale(x_a)), by = uid]
envelopes[, x_b := as.numeric(scale(x_b)), by = uid]
series <- split(envelopes[, .(uid, x_a, x_b)], by = "uid", keep.by = FALSE)
rm(envelopes)
fs <- round(median(windows$n_samples / windows$duration_s), 3)
mindiagline <- max(2L, as.integer(round(line_s * fs)))
cat(paste(sessions, collapse = ", "), ":", nrow(windows), "windows at", fs, "Hz,", cores, "cores\n")
t0 <- Sys.time()
elapsed <- function() sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

# --- Radius per area: median over windows of the target quantile of |A - B| -----------------------------------
q <- unlist(mclapply(windows$uid, function(id) {
  w <- series[[id]]
  quantile(abs(outer(w$x_a, w$x_b, "-")), target / 100, names = FALSE)
}, mc.cores = cores))
parameters <- windows[, .(radius = median(q[.I]), n_windows = .N), by = area]
parameters[, `:=`(embed = embed, delay = delay, target = target, fs = fs, mindiagline = mindiagline,
                  mindiagline_s = mindiagline / fs, sessions = paste(sessions, collapse = " "))]
print(parameters)
fwrite(parameters, file.path(out, "parameters.csv"))

# --- cRQA of every window and of the surrogate pairs ----------------------------------------------------------
band <- round(band_s * fs)
max_prof <- round(profile_s * fs)
measure <- function(xa, xb, radius) {
  r <- crqa(xa, xb, delay = delay, embed = embed, radius = radius, rescale = 0, normalize = 0,
            mindiagline = mindiagline, minvertline = mindiagline, tw = 0, whiteline = FALSE, recpt = FALSE,
            side = "both", method = "crqa", metric = "euclidean", datatype = "continuous")
  n <- nrow(r$RP)
  nz <- summary(r$RP)                       # recurrent points: i = A's time, j = B's time
  lag <- nz$j - nz$i                        # positive = B's state matches A's earlier state (B follows A)
  counts <- tabulate(lag[abs(lag) <= max_prof] + max_prof + 1, nbins = 2 * max_prof + 1)
  cells <- n - abs(-max_prof:max_prof)      # cells on each diagonal
  in_band <- abs(-max_prof:max_prof) <= band
  list(values = data.table(RR = r$RR, DET = r$DET, L_s = r$L / fs, LAM = r$LAM,
                           lagRR_2s = 100 * sum(counts[in_band]) / sum(cells[in_band])),
       profile = 100 * counts / cells)
}

joint <- windows[condition == "joint", .(uid, session, area, pair, round, trial)]
surrogates <- merge(joint, joint, by = c("session", "area", "pair", "trial"), suffixes = c("_a", "_b"),
                    allow.cartesian = TRUE)[round_a != round_b]
jobs <- rbind(windows[, .(kind = "real", uid, uid_b = uid, area)],
              surrogates[, .(kind = "surrogate", uid = uid_a, uid_b, area)])
radius_of <- setNames(parameters$radius, parameters$area)
cat(nrow(jobs), "cRQA jobs (", sum(jobs$kind == "real"), "real,", sum(jobs$kind == "surrogate"), "surrogate ) |",
    elapsed(), "\n")

computed <- mclapply(seq_len(nrow(jobs)), function(k) {
  j <- jobs[k]
  xa <- series[[j$uid]]$x_a
  xb <- series[[j$uid_b]]$x_b
  n <- min(length(xa), length(xb))
  m <- measure(xa[1:n], xb[1:n], radius_of[[j$area]])
  list(values = cbind(j, n_samples = n, m$values),
       profile = if (j$kind == "real") data.table(uid = j$uid, lag_s = (-max_prof:max_prof) / fs, RR = m$profile))
}, mc.cores = cores)
failed <- sapply(computed, inherits, "try-error")
if (any(failed)) stop(sum(failed), " cRQA jobs failed, e.g. ", computed[[which(failed)[1]]])
results <- rbindlist(lapply(computed, `[[`, "values"))
results <- merge(results, windows[, .(uid, window_id, session, pair, participant_a, participant_b, condition, round,
                                      trial, role, duration_s, bad_fraction_a, bad_fraction_b)], by = "uid")
results[kind == "surrogate", role := "surrogate"]
fwrite(windows, file.path(out, "windows.csv"))
fwrite(results, file.path(out, "results.csv"))
fwrite(rbindlist(lapply(computed, `[[`, "profile")), file.path(out, "profiles.csv"))

print(results[, .(n = .N, RR = median(RR), DET = median(DET), lagRR_2s = median(lagRR_2s)),
              by = .(area, role)][order(area, role)], digits = 3)
cat("saved to", out, "|", elapsed(), "\n")
