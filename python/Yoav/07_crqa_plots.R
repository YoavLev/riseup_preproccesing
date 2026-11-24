# 07_crqa_plots.R: figures of the cRQA results of one session folder, after 06_crqa.R. For viewing and discussing
# the method (a demo of one pair so far), not a statistical test.
#
# Figures (in derivatives/crqa/<folder>/plots/), per area (visual, motor):
#   1. series_<area>.png   the two z-scored envelopes of one window per role (first 60 s): what cRQA compares.
#   2. crp_<area>.png      cross-recurrence plots of the same windows (first 60 s). A dot at (t_A, t_B): the state of
#                          A at t_A is close to the state of B at t_B (within the radius). Diagonal lines = A and B
#                          go through the same sequence of states; lines on the main diagonal = at the same time.
#   3. profile_<area>.png  diagonal cross-recurrence profile: recurrence rate as a function of the lag between A and
#                          B (-5..+5 s), one thin line per window, thick = mean per role. A peak at lag 0 = the two
#                          envelopes are similar at the same moment; a peak at +x s = B follows A by x s.
#   4. measures.png        %REC (RR) and %DET of every window, by role and area (as saved by 06_crqa.R).
#
# The cRQA parameters are those 06_crqa.R used (crqa_parameters_<folder>.csv, used = TRUE), the envelopes are
# z-scored per window as in 06_crqa.R, so the plots show exactly what the measures were computed on.
#
# Usage (on UCloud, after scripts/install_r.sh):
#   Rscript 07_crqa_plots.R session03 [cores]
#   Rscript 07_crqa_plots.R session03_smooth1hz [cores]

suppressMessages({
  library(crqa)
  library(data.table)
  library(parallel)
  library(Matrix)
})
# data.table's own threads (OpenMP) make the forked workers of mclapply() hang (they wait forever on a lock):
# one thread, the parallel work is done by mclapply() alone
setDTthreads(1)

args <- commandArgs(trailingOnly = TRUE)
folder <- if (length(args) >= 1) args[1] else "session03"
cores <- if (length(args) >= 2) as.integer(args[2]) else 8
dir <- file.path(Sys.getenv("RT_CRQA_ROOT", "/work/YoavLevy#3923/derivatives/crqa"), folder)
out <- file.path(dir, "plots")
dir.create(out, showWarnings = FALSE)

# --- Settings -------------------------------------------------------------------------------------------------
show_s <- 60           # seconds of each example window shown in figures 1 and 2
max_lag_s <- 5         # lag range of the diagonal profile [s], both sides
roles <- c("baseline", "solo", "teammates", "opponents")
role_colour <- c(baseline = "#999999", solo = "#E69F00", teammates = "#0072B2", opponents = "#D55E00")

# --- Data -----------------------------------------------------------------------------------------------------
windows <- fread(file.path(dir, "windows.csv"))
results <- fread(file.path(dir, paste0("crqa_", folder, ".csv")))
params <- fread(file.path(dir, paste0("crqa_parameters_", folder, ".csv")))[used == TRUE]
envelopes <- fread(file.path(dir, "envelopes.csv"))
envelopes[, x_a := as.numeric(scale(x_a)), by = window_id]   # z-score per window and person, as in 06_crqa.R
envelopes[, x_b := as.numeric(scale(x_b)), by = window_id]
series <- split(envelopes, by = "window_id")
fs <- round(median(windows$n_samples / windows$duration_s))   # envelope sampling rate (25 Hz)
max_lag <- max_lag_s * fs
windows[, role := factor(role, levels = roles)]
results[, role := factor(role, levels = roles)]

# One example window per area and role: the first trial of that role (lowest round, then trial)
examples <- windows[order(area, role, round, trial), .SD[1], by = .(area, role)]

# --- cRQA of every window: the recurrence matrix, for the profile (and the plot of the examples) ----------------
run_crqa <- function(id) {
  w <- series[[as.character(id)]]
  p <- params[area == windows[window_id == id, area]]
  crqa(w$x_a, w$x_b, delay = p$delay, embed = p$embed, radius = p$radius, rescale = 0, normalize = 0,
       mindiagline = 2, minvertline = 2, tw = 0, whiteline = FALSE, recpt = FALSE, side = "both",
       method = "crqa", metric = "euclidean", datatype = "continuous")
}

# Diagonal profile from the recurrence matrix: for lag k, the share of recurrent points on the diagonal
# t_B - t_A = k (in %). Positive lag = B's state matches A's state k samples later (B follows A).
diagonal_profile <- function(rp) {
  n <- nrow(rp)
  nz <- summary(rp)                        # row (i = A's time), column (j = B's time) of every recurrent point
  lag <- nz$j - nz$i
  lag <- lag[abs(lag) <= max_lag]
  counts <- tabulate(lag + max_lag + 1, nbins = 2 * max_lag + 1)
  100 * counts / (n - abs(-max_lag:max_lag))
}

cat(folder, ":", nrow(windows), "windows, cRQA for the profiles on", cores, "cores\n")
computed <- mclapply(windows$window_id, function(id) {
  r <- run_crqa(id)
  keep <- id %in% examples$window_id
  list(window_id = id, RR = r$RR, profile = diagonal_profile(r$RP),
       rp = if (keep) r$RP[seq_len(min(nrow(r$RP), show_s * fs)), seq_len(min(ncol(r$RP), show_s * fs))] else NULL)
}, mc.cores = cores)
names(computed) <- sapply(computed, `[[`, "window_id")

# Check: the recurrence rate recomputed here must match what 06_crqa.R saved
check <- merge(results[, .(window_id, RR_06 = RR)], data.table(window_id = as.integer(names(computed)),
               RR_07 = sapply(computed, `[[`, "RR")), by = "window_id")
cat("largest RR difference to 06_crqa.R:", signif(max(abs(check$RR_06 - check$RR_07)), 3), "%\n")

title_of <- function(w) {
  r <- results[window_id == w$window_id]
  sprintf("%s: %s r%d t%d   (whole window %.0f s: RR %.1f %%, DET %.1f %%)", as.character(w$role), w$condition, w$round,
          w$trial, w$duration_s, r$RR, r$DET)
}

for (ar in unique(windows$area)) {
  ex <- examples[area == ar]
  p <- params[area == ar]
  pair <- ex$pair[1]
  setting <- sprintf("%s, %s %s | delay %d (%.2f s), embedding %d, radius %.2f", pair, ar, ex$band[1], p$delay,
                     p$delay / fs, p$embed, p$radius)

  # --- 1. the envelopes -------------------------------------------------------------------------------------
  png(file.path(out, paste0("series_", ar, ".png")), width = 1400, height = 250 * nrow(ex) + 80, res = 110)
  par(mfrow = c(nrow(ex), 1), mar = c(3, 4, 2, 1), oma = c(0, 0, 2.5, 0), mgp = c(2, 0.6, 0))
  for (k in seq_len(nrow(ex))) {
    w <- series[[as.character(ex$window_id[k])]][sample < show_s * fs]
    t <- w$sample / fs
    plot(t, w$x_a, type = "l", col = "black", lwd = 1, ylim = range(c(w$x_a, w$x_b)),
         xlab = "time in the window (s)", ylab = "envelope (z)", main = title_of(ex[k]), cex.main = 0.95)
    lines(t, w$x_b, col = role_colour[as.character(ex$role[k])], lwd = 1)
    legend("topright", c(ex$participant_a[k], ex$participant_b[k]), col = c("black", role_colour[as.character(ex$role[k])]),
           lwd = 1, bty = "n", cex = 0.8, horiz = TRUE)
  }
  mtext(paste("Envelopes compared by cRQA (first", show_s, "s;", setting, ")"), outer = TRUE, cex = 0.9)
  dev.off()

  # --- 2. cross-recurrence plots ----------------------------------------------------------------------------
  png(file.path(out, paste0("crp_", ar, ".png")), width = 1400, height = 1450, res = 110)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2.5, 0), mgp = c(2.2, 0.6, 0), pty = "s")
  for (k in seq_len(nrow(ex))) {
    rp <- as.matrix(computed[[as.character(ex$window_id[k])]]$rp)
    t <- (seq_len(nrow(rp)) - 1) / fs
    image(t, t, rp * 1, col = c("white", role_colour[as.character(ex$role[k])]), useRaster = TRUE,
          xlab = paste("time of", ex$participant_a[k], "(s)"), ylab = paste("time of", ex$participant_b[k], "(s)"),
          main = sub("   ", "\n", title_of(ex[k])), cex.main = 0.9)
    abline(0, 1, col = "grey40", lty = 3)
  }
  mtext(paste("Cross-recurrence plots (first", show_s, "s;", setting, ")"), outer = TRUE, cex = 0.9)
  dev.off()

  # --- 3. diagonal profile ----------------------------------------------------------------------------------
  lags <- (-max_lag:max_lag) / fs
  ids <- windows[area == ar, window_id]
  prof <- do.call(rbind, lapply(as.character(ids), function(id) computed[[id]]$profile))
  role_of <- windows[match(ids, window_id), role]
  png(file.path(out, paste0("profile_", ar, ".png")), width = 1300, height = 800, res = 110)
  par(mar = c(4, 4, 3, 1), mgp = c(2.2, 0.6, 0))
  plot(NA, xlim = range(lags), ylim = range(prof), xlab = "lag (s): positive = B follows A",
       ylab = "recurrence rate on the diagonal (%)",
       main = paste("Diagonal cross-recurrence profile\n", setting), cex.main = 0.9)
  abline(v = 0, col = "grey80")
  for (i in seq_len(nrow(prof))) lines(lags, prof[i, ], col = adjustcolor(role_colour[as.character(role_of[i])], 0.25))
  n_role <- table(droplevels(role_of))
  for (ro in names(n_role)) lines(lags, colMeans(prof[role_of == ro, , drop = FALSE]), col = role_colour[ro], lwd = 3)
  legend("topright", sprintf("%s (n = %d)", names(n_role), n_role), col = role_colour[names(n_role)], lwd = 3, bty = "n")
  dev.off()
  fwrite(data.table(window_id = rep(ids, each = length(lags)), role = rep(role_of, each = length(lags)),
                    lag_s = rep(lags, length(ids)), RR = as.vector(t(prof))),
         file.path(out, paste0("profile_", ar, ".csv")))
}

# --- 4. %REC and %DET by role --------------------------------------------------------------------------------
png(file.path(out, "measures.png"), width = 1400, height = 1000, res = 110)
areas <- unique(results$area)
par(mfrow = c(2, length(areas)), mar = c(4, 4, 3, 1), mgp = c(2.2, 0.6, 0))
set.seed(97)   # the horizontal jitter of the points
for (measure in c("RR", "DET")) for (ar in areas) {
  d <- results[area == ar]
  x <- as.integer(d$role)
  plot(x + runif(nrow(d), -0.15, 0.15), d[[measure]], pch = 19, col = role_colour[as.character(d$role)],
       xlim = c(0.5, length(roles) + 0.5), xaxt = "n", xlab = "", ylab = paste0(measure, " (%)"),
       main = sprintf("%s: %s (%s)", ar, if (measure == "RR") "%REC" else "%DET", d$pair[1]))
  axis(1, at = seq_along(roles), labels = roles)
  med <- d[, .(m = median(get(measure))), by = role]
  segments(as.integer(med$role) - 0.3, med$m, as.integer(med$role) + 0.3, med$m, lwd = 3)
}
dev.off()

print(results[, .(n = .N, RR = median(RR), DET = median(DET), L = median(L)), by = .(area, role)][order(area, role)],
      digits = 3)
cat("figures saved to", out, "\n")
