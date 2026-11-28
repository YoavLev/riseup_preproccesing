# 07_crqa_plots.R: figures and tables of the final cRQA (06_crqa.R), all sessions together.
#
# Within-pair differences: for each pair, the mean over its windows of one role minus the mean over its windows of
# the other role (only pairs that had both). Contrasts:
#   teammates - opponents   (RQ1; both roles occur in every joint block, so not confounded with time)
#   teammates - solo        (RQ2; Solo always comes first: confounded with time / order)
#   joint - surrogate       (real joint pairs minus the same people in different blocks: coupling of the moment)
# Pairs share participants, so they are not independent: the CIs (t over pairs) and Wilcoxon p are a guide.
# Mixed model (if lme4 is installed), joint windows only:
#   measure ~ role + block + (1|session) + (1|pair) + (1|participant_a) + (1|participant_b)
#   role effect tested with a likelihood-ratio test against the model without role (ML fit).
#
# Usage: Rscript 07_crqa_plots.R [folder]   (default final_smooth1hz_5hz)
# Output: derivatives/crqa/<folder>/summary/  contrasts.png, contrasts_RR.png, roles.png, profile.png,
#         sessions.png, crp_example_<area>.png, contrasts.csv, contrasts_per_pair.csv, models.csv, counts.csv,
#         outliers.csv

suppressMessages({
  library(crqa)
  library(data.table)
  library(Matrix)
})
setDTthreads(1)
args <- commandArgs(trailingOnly = TRUE)
folder <- if (length(args) >= 1) args[1] else "final_smooth1hz_5hz"
root <- Sys.getenv("RT_CRQA_ROOT", "/work/YoavLevy#3923/derivatives/crqa")
dir <- file.path(root, folder)
out <- file.path(dir, "summary")
dir.create(out, showWarnings = FALSE)

roles <- c("baseline", "solo", "teammates", "opponents", "surrogate")
role_colour <- c(baseline = "#999999", solo = "#E69F00", teammates = "#0072B2", opponents = "#D55E00",
                 surrogate = "#555555")
results <- fread(file.path(dir, "results.csv"))
params <- fread(file.path(dir, "parameters.csv"))
sessions <- sort(unique(results$session))
session_colour <- setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#999999",
                             "#000000")[seq_along(sessions)], sessions)
measures <- c("lagRR_2s", "RR", "DET")
areas <- sort(unique(results$area))   # "visual" (only occipital alpha is analysed)

# --- Counts ---------------------------------------------------------------------------------------------------
counts <- results[, .(windows = .N, pairs = uniqueN(paste(session, pair))), by = .(area, role)]
fwrite(counts, file.path(out, "counts.csv"))

# --- Within-pair contrasts ------------------------------------------------------------------------------------
pm <- results[, lapply(.SD, mean), by = .(session, area, pair, role), .SDcols = measures]
pm <- rbind(pm, results[kind == "real" & condition == "joint", lapply(.SD, mean), by = .(session, area, pair),
                        .SDcols = measures][, role := "joint"])
contrast <- function(r1, r2) {
  w <- dcast(pm[role %in% c(r1, r2)], session + area + pair ~ role, value.var = measures)
  rbindlist(lapply(measures, function(m) {
    w[, .(measure = m, contrast = paste(r1, "-", r2), diff = get(paste0(m, "_", r1)) - get(paste0(m, "_", r2))),
      by = .(session, area, pair)]
  }))[!is.na(diff)]
}
diffs <- rbind(contrast("teammates", "opponents"), contrast("teammates", "solo"), contrast("joint", "surrogate"))
summ <- diffs[, {
  n <- .N; m <- mean(diff); se <- sd(diff) / sqrt(n)
  .(n_pairs = n, mean = m, lo = m - qt(0.975, n - 1) * se, hi = m + qt(0.975, n - 1) * se,
    share_positive = mean(diff > 0), wilcoxon_p = if (n >= 5) suppressWarnings(wilcox.test(diff)$p.value) else NA)
}, by = .(contrast, measure, area)]
fwrite(summ, file.path(out, "contrasts.csv"))
fwrite(diffs, file.path(out, "contrasts_per_pair.csv"))
print(summ, digits = 3)

# --- Mixed models (joint windows) -----------------------------------------------------------------------------
if (requireNamespace("lme4", quietly = TRUE)) {
  d <- results[kind == "real" & condition == "joint"]
  d[, `:=`(role = factor(role, levels = c("opponents", "teammates")), block = factor(round),
           pair_id = paste(session, pair))]
  models <- rbindlist(lapply(areas, function(ar) rbindlist(lapply(c("lagRR_2s", "RR"), function(m) {
    dd <- d[area == ar]
    f1 <- as.formula(paste(m, "~ role + block + (1|session) + (1|pair_id) + (1|participant_a) + (1|participant_b)"))
    f0 <- update(f1, . ~ . - role)
    m1 <- suppressMessages(lme4::lmer(f1, data = dd, REML = FALSE))
    m0 <- suppressMessages(lme4::lmer(f0, data = dd, REML = FALSE))
    lrt <- anova(m0, m1)
    co <- summary(m1)$coefficients["roleteammates", ]
    data.table(area = ar, measure = m, n_windows = nrow(dd), teammates_minus_opponents = co[["Estimate"]],
               se = co[["Std. Error"]], chisq = lrt$Chisq[2], p_lrt = lrt$`Pr(>Chisq)`[2],
               singular = lme4::isSingular(m1))
  }))))
  fwrite(models, file.path(out, "models.csv"))
  print(models, digits = 3)
} else cat("lme4 not installed: no mixed models\n")

# --- Figure: contrasts (one figure per measure) ---------------------------------------------------------------
contrast_titles <- c("teammates - opponents" = "Teammates minus opponents",
                     "teammates - solo" = "Teammates minus playing alone",
                     "joint - surrogate" = "Real pairs minus shuffled pairs")
plot_contrasts <- function(m, file, ylab) {
  png(file.path(out, file), width = 1800, height = 1100, res = 130)
  par(mfrow = c(length(areas), 3), mar = c(3, 4.5, 3, 1), mgp = c(2.6, 0.7, 0))
  set.seed(97)
  for (ar in areas) for (ct in names(contrast_titles)) {
    d <- diffs[measure == m & area == ar & contrast == ct]
    s <- summ[measure == m & area == ar & contrast == ct]
    lim <- range(c(d$diff, s$lo, s$hi, 0))
    plot(runif(nrow(d), 0.75, 1.25), d$diff, pch = 19, cex = 0.9, col = adjustcolor(session_colour[d$session], 0.7),
         xlim = c(0.4, 1.9), ylim = lim, xaxt = "n", xlab = "", ylab = ylab,
         main = sprintf("%s: %s", ar, contrast_titles[[ct]]), cex.main = 1)
    abline(h = 0, lty = 2, col = "grey50")
    sig <- s$lo > 0 | s$hi < 0
    arrows(1.6, s$lo, 1.6, s$hi, angle = 90, code = 3, length = 0.05, lwd = 2, col = if (sig) "#0072B2" else "black")
    points(1.6, s$mean, pch = 19, cex = 1.6, col = if (sig) "#0072B2" else "black")
    mtext(sprintf("%d pairs, mean %.2f", s$n_pairs, s$mean), side = 1, line = 1, cex = 0.7)
    if (ar == areas[1] && ct == names(contrast_titles)[1])
      legend("topleft", sessions, col = session_colour, pch = 19, bty = "n", cex = 0.75)
  }
  dev.off()
}
plot_contrasts("lagRR_2s", "contrasts.png", "difference in recurrence within 2 s (points)")
plot_contrasts("RR", "contrasts_RR.png", "difference in whole-plot recurrence (points)")

# --- Figure: lagRR_2s by role ---------------------------------------------------------------------------------
png(file.path(out, "roles.png"), width = 1700, height = 700, res = 130)
par(mfrow = c(1, length(areas)), mar = c(3, 4.5, 3, 1), mgp = c(2.6, 0.7, 0))
for (ar in areas) {
  d <- pm[area == ar & role %in% roles]
  d[, x := match(role, roles)]
  plot(d$x + runif(nrow(d), -0.2, 0.2), d$lagRR_2s, pch = 19, cex = 0.7,
       col = adjustcolor(role_colour[d$role], 0.6), xlim = c(0.5, length(roles) + 0.5), xaxt = "n", xlab = "",
       ylab = "recurrence within 2 s (%), mean per pair", main = ar)
  axis(1, at = seq_along(roles), labels = c("bubbles", "alone", "teammates", "opponents", "shuffled"))
  med <- d[, .(v = median(lagRR_2s)), by = x]
  segments(med$x - 0.3, med$v, med$x + 0.3, med$v, lwd = 3)
}
dev.off()

# --- Figure: per session, teammates - opponents ---------------------------------------------------------------
png(file.path(out, "sessions.png"), width = 1500, height = 650, res = 130)
par(mfrow = c(1, length(areas)), mar = c(3, 4.5, 3, 1), mgp = c(2.6, 0.7, 0))
for (ar in areas) {
  d <- diffs[measure == "lagRR_2s" & area == ar & contrast == "teammates - opponents"]
  s <- d[, .(m = mean(diff), se = sd(diff) / sqrt(.N), n = .N), by = session][order(session)]
  x <- seq_len(nrow(s))
  plot(x, s$m, pch = 19, col = session_colour[s$session], cex = 1.5, xlim = c(0.5, nrow(s) + 0.5),
       ylim = range(c(s$m - 2 * s$se, s$m + 2 * s$se, 0), na.rm = TRUE), xaxt = "n", xlab = "",
       ylab = "teammates - opponents (points)", main = sprintf("%s: per session (mean +- 2 SE)", ar))
  arrows(x, s$m - 2 * s$se, x, s$m + 2 * s$se, angle = 90, code = 3, length = 0.04)
  abline(h = 0, lty = 2, col = "grey50")
  axis(1, at = x, labels = sprintf("%s\n(%d pairs)", sub("session", "S", s$session), s$n), padj = 0.5, cex.axis = 0.85)
}
dev.off()

# --- Figure: diagonal profile by role -------------------------------------------------------------------------
prof <- merge(fread(file.path(dir, "profiles.csv")), fread(file.path(dir, "windows.csv"))[, .(uid, role, area)],
              by = "uid")
png(file.path(out, "profile.png"), width = 1700, height = 700, res = 130)
par(mfrow = c(1, length(areas)), mar = c(4, 4.5, 3, 1), mgp = c(2.6, 0.7, 0))
for (ar in areas) {
  p <- prof[area == ar, .(RR = mean(RR), se = sd(RR) / sqrt(.N)), by = .(role, lag_s)]
  p[is.na(se), se := 0]
  plot(NA, xlim = range(p$lag_s), ylim = range(c(p$RR - 2 * p$se, p$RR + 2 * p$se)),
       xlab = "time shift between the two people (s)", ylab = "recurrence (%)", main = ar)
  rect(-2, -1e3, 2, 1e3, col = "grey93", border = NA)
  for (ro in intersect(roles, unique(p$role))) {
    q <- p[role == ro][order(lag_s)]
    polygon(c(q$lag_s, rev(q$lag_s)), c(q$RR - 2 * q$se, rev(q$RR + 2 * q$se)),
            col = adjustcolor(role_colour[ro], 0.15), border = NA)
    lines(q$lag_s, q$RR, col = role_colour[ro], lwd = 2)
  }
  legend("topright", c("bubbles", "alone", "teammates", "opponents")[match(intersect(roles, unique(p$role)), roles)],
         col = role_colour[intersect(roles, unique(p$role))], lwd = 2, bty = "n", cex = 0.85)
}
dev.off()

# --- Figure: a cross-recurrence plot (first 60 s of one teammate window) --------------------------------------
w <- fread(file.path(dir, "windows.csv"))
for (ar in areas) {
  ex <- w[area == ar & role == "teammates"][order(session, pair, round, trial)][1]
  e <- fread(file.path(root, paste0(ex$session, "_", sub("^final_", "", folder)), "envelopes.csv"))[window_id ==
                                                                                                         ex$window_id]
  xa <- as.numeric(scale(e$x_a)); xb <- as.numeric(scale(e$x_b))
  p <- params[area == ar]
  r <- crqa(xa, xb, delay = 1, embed = 1, radius = p$radius, rescale = 0, normalize = 0, mindiagline = p$mindiagline,
            minvertline = p$mindiagline, tw = 0, whiteline = FALSE, recpt = FALSE, side = "both", method = "crqa",
            metric = "euclidean", datatype = "continuous")
  n <- min(nrow(r$RP), round(60 * p$fs))
  t <- (seq_len(n) - 1) / p$fs
  png(file.path(out, sprintf("crp_example_%s.png", ar)), width = 1500, height = 800, res = 130)
  layout(matrix(c(1, 2), 1), widths = c(1.3, 1))
  par(mar = c(4, 4.5, 3, 1), mgp = c(2.6, 0.7, 0))
  plot(t, xa[1:n], type = "l", col = "#0072B2", lwd = 2, xlab = "time (s)", ylab = "envelope (z)",
       ylim = range(c(xa[1:n], xb[1:n])), main = "two teammates' envelopes (60 s)")
  lines(t, xb[1:n], col = "#D55E00", lwd = 2)
  par(pty = "s")
  image(t, t, as.matrix(r$RP[1:n, 1:n]) * 1, col = c("white", "black"), useRaster = TRUE,
        xlab = "time, person A (s)", ylab = "time, person B (s)", main = "recurrence plot")
  abline(0, 1, col = "grey60", lty = 3)
  abline(2, 1, col = "#0072B2", lty = 2); abline(-2, 1, col = "#0072B2", lty = 2)
  dev.off()
}

# --- Outlier windows (flagged, not removed): > 3 IQR from the median of their area and role, or > 20 % bad ------
d <- results[kind == "real"]
d[, outlier := abs(lagRR_2s - median(lagRR_2s)) > 3 * IQR(lagRR_2s) | abs(RR - median(RR)) > 3 * IQR(RR),
  by = .(area, role)]
d[, bad := pmax(bad_fraction_a, bad_fraction_b)]
fwrite(d[outlier | bad > 0.2, .(session, area, role, pair, condition, round, trial, duration_s, RR, lagRR_2s,
                                bad_fraction_a, bad_fraction_b, outlier)], file.path(out, "outliers.csv"))
cat("outlier/bad windows:", d[outlier | bad > 0.2, .N], "of", nrow(d), "\n")
cat("saved to", out, "\n")
