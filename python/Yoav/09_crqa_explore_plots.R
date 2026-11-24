# 09_crqa_explore_plots.R: figures and summary tables of the parameter exploration (08_crqa_explore.R), over all
# window folders of the given sessions (<session>_smooth<X>hz_<Y>hz). For viewing and discussion, not the final test.
#
# Figures (derivatives/crqa/explore_summary/):
#   curves.png           pooled AMI and FNN curves per smoothing level and area, with the chosen delay
#   roles_<setting>.png  RR, DET and lagRR_2s per window, by role (and surrogate pairs), per area
#   contrasts.png        per setting: within-pair differences (teammates - opponents, teammates - solo, real joint -
#                        surrogate) of each measure, mean over pairs with 95 % CI. The robustness overview.
#   profile_<setting>.png diagonal profile (RR per lag) by role, mean over windows
#   crp_<setting>_<area>.png cross-recurrence plots (first 60 s) of one pair, one window per role
# Tables: parameters_all.csv, contrasts.csv, outliers.csv
#
# Within-pair differences: for each pair, the mean over its windows of one role minus the mean over its windows
# of the other role (pairs that had both roles). Pairs share participants, so they are not independent: the CIs
# (t over pairs) are only a guide.
#
# Usage: Rscript 09_crqa_explore_plots.R session03,session04 smooth1hz_5hz/fnn/5 [cores]
#   2nd argument = the setting shown in the role, profile and CRP figures (folder suffix / embedding / target)

suppressMessages({
  library(crqa)
  library(data.table)
  library(Matrix)
})
setDTthreads(1)

args <- commandArgs(trailingOnly = TRUE)
sessions <- strsplit(if (length(args) >= 1) args[1] else "session03,session04", ",")[[1]]
main <- strsplit(if (length(args) >= 2) args[2] else "smooth1hz_5hz/fnn/5", "/")[[1]]
root <- Sys.getenv("RT_CRQA_ROOT", "/work/YoavLevy#3923/derivatives/crqa")
out <- file.path(root, "explore_summary")
dir.create(out, showWarnings = FALSE)
roles <- c("baseline", "solo", "teammates", "opponents", "surrogate")
role_colour <- c(baseline = "#999999", solo = "#E69F00", teammates = "#0072B2", opponents = "#D55E00",
                 surrogate = "#555555")

# --- Read every folder --------------------------------------------------------------------------------------
folders <- unlist(lapply(sessions, function(s) list.files(root, pattern = paste0("^", s, "_smooth.*hz_.*hz$"))))
folders <- folders[file.exists(file.path(root, folders, "explore", "results.csv"))]
cat("folders:", paste(folders, collapse = ", "), "\n")
read_all <- function(file) rbindlist(lapply(folders, function(f) {
  d <- fread(file.path(root, f, "explore", file))
  d[, `:=`(folder = f, smoothing = sub("^session\\d+_", "", f))]
}), fill = TRUE)
params <- read_all("parameters.csv")
curves <- read_all("curves.csv")
results <- read_all("results.csv")
results[, setting := paste(smoothing, embedding, target, sep = "/")]
fwrite(params, file.path(out, "parameters_all.csv"))
smooth_levels <- unique(params[order(fs), smoothing])

# --- 1. AMI and FNN curves ----------------------------------------------------------------------------------
png(file.path(out, "curves.png"), width = 1500, height = 420 * length(smooth_levels), res = 110)
par(mfrow = c(length(smooth_levels), 2), mar = c(4, 4, 3, 1), mgp = c(2.2, 0.6, 0))
area_col <- c(visual = "#009E73", motor = "#CC79A7")
for (sm in smooth_levels) {
  for (kind in c("AMI", "FNN")) {
    d <- curves[smoothing == sm & curve == kind]
    xvar <- if (kind == "AMI") "x_s" else "x"
    plot(NA, xlim = range(d[[xvar]]), ylim = range(d$y, na.rm = TRUE), xlab = if (kind == "AMI") "lag (s)" else
         "embedding dimension", ylab = if (kind == "AMI") "mutual information" else "false nearest neighbours (%)",
         main = sprintf("%s: %s (average over windows and people)", sm, kind), cex.main = 0.95)
    for (f in unique(d$folder)) for (ar in unique(d$area)) {
      dd <- d[folder == f & area == ar]
      lines(dd[[xvar]], dd$y, col = area_col[ar], lty = if (grepl("session03", f)) 1 else 2, lwd = 1.5)
      p <- params[folder == f & area == ar][1]
      if (kind == "AMI") points(p$delay_s, dd[x == p$delay, y], pch = 19, col = area_col[ar])
    }
    legend("topright", c(paste(names(area_col)), "session03", "session04"), col = c(area_col, "black", "black"),
           lty = c(1, 1, 1, 2), bty = "n", cex = 0.8)
  }
}
dev.off()

# --- 2. Within-pair contrasts for every setting --------------------------------------------------------------
measures <- c("RR", "DET", "lagRR_2s")
pair_means <- results[, lapply(.SD, mean), by = .(session, setting, smoothing, embedding, target, area, pair, role),
                      .SDcols = measures]
pair_means[role %in% c("teammates", "opponents"), joint := TRUE]
joint_real <- results[kind == "real" & condition == "joint",
                      lapply(.SD, mean), by = .(session, setting, smoothing, embedding, target, area, pair),
                      .SDcols = measures][, role := "joint"]
pair_means <- rbind(pair_means[, !"joint"], joint_real)
contrast <- function(r1, r2) {
  w <- dcast(pair_means[role %in% c(r1, r2)], session + setting + smoothing + embedding + target + area + pair ~ role,
             value.var = measures)
  rbindlist(lapply(measures, function(m) {
    w[, .(measure = m, contrast = paste(r1, "-", r2), diff = get(paste0(m, "_", r1)) - get(paste0(m, "_", r2))),
      by = .(session, setting, smoothing, embedding, target, area, pair)]
  }))[!is.na(diff)]
}
diffs <- rbind(contrast("teammates", "opponents"), contrast("teammates", "solo"), contrast("joint", "surrogate"))
summ <- diffs[, {
  n <- .N; m <- mean(diff); se <- sd(diff) / sqrt(n)
  p <- if (n >= 5) suppressWarnings(wilcox.test(diff)$p.value) else NA_real_
  .(n_pairs = n, mean = m, lo = m - qt(0.975, n - 1) * se, hi = m + qt(0.975, n - 1) * se,
    positive = mean(diff > 0), wilcoxon_p = p)
}, by = .(contrast, measure, area, smoothing, embedding, target, setting)]
fwrite(summ, file.path(out, "contrasts.csv"))
fwrite(diffs, file.path(out, "contrasts_per_pair.csv"))

settings <- unique(summ[order(match(smoothing, smooth_levels), embedding, target), setting])
contrasts <- unique(summ$contrast)
png(file.path(out, "contrasts.png"), width = 1700, height = 330 * length(contrasts) * 2, res = 110)
par(mfrow = c(length(contrasts) * 2, length(measures)), mar = c(7.5, 4, 2.5, 1), mgp = c(2.2, 0.6, 0))
for (ct in contrasts) for (ar in c("visual", "motor")) for (m in measures) {
  d <- summ[contrast == ct & area == ar & measure == m][match(settings, setting)]
  x <- seq_along(settings)
  plot(x, d$mean, ylim = range(c(d$lo, d$hi, 0), na.rm = TRUE), pch = 19, xaxt = "n", xlab = "",
       ylab = paste("difference in", m, "(points)"), main = sprintf("%s, %s: %s", ar, m, ct), cex.main = 0.95,
       col = ifelse(d$lo > 0 | d$hi < 0, "#0072B2", "black"))
  arrows(x, d$lo, x, d$hi, angle = 90, code = 3, length = 0.03)
  abline(h = 0, col = "grey60", lty = 2)
  axis(1, at = x, labels = settings, las = 2, cex.axis = 0.7)
}
dev.off()

# --- 3. Main setting: measures by role ----------------------------------------------------------------------
main_setting <- paste(main, collapse = "/")
tag <- gsub("/", "_", main_setting)
d <- results[setting == main_setting]
d[, role := factor(role, levels = roles)]
png(file.path(out, paste0("roles_", tag, ".png")), width = 1600, height = 950, res = 110)
par(mfrow = c(2, length(measures)), mar = c(4, 4, 3, 1), mgp = c(2.2, 0.6, 0))
set.seed(97)
for (ar in c("visual", "motor")) for (m in measures) {
  dd <- d[area == ar]
  x <- as.integer(dd$role)
  plot(x + runif(nrow(dd), -0.18, 0.18), dd[[m]], pch = ifelse(dd$session == "session03", 19, 17), cex = 0.6,
       col = adjustcolor(role_colour[as.character(dd$role)], 0.6), xlim = c(0.5, length(roles) + 0.5), xaxt = "n",
       xlab = "", ylab = paste0(m, " (%)"), main = sprintf("%s: %s  (%s)", ar, m, main_setting), cex.main = 0.95)
  axis(1, at = seq_along(roles), labels = roles, cex.axis = 0.85)
  med <- dd[, .(v = median(get(m))), by = role]
  segments(as.integer(med$role) - 0.3, med$v, as.integer(med$role) + 0.3, med$v, lwd = 3)
}
dev.off()

# --- 4. Main setting: diagonal profile by role --------------------------------------------------------------
prof <- rbindlist(lapply(folders[grepl(paste0("_", main[1], "$"), folders)], function(f) {
  p <- fread(file.path(root, f, "explore", "profiles.csv"))[embedding == main[2] & target == as.numeric(main[3])]
  w <- fread(file.path(root, f, "windows.csv"))[, .(window_id, role, area)]
  merge(p, w, by = "window_id")[, folder := f]
}))
png(file.path(out, paste0("profile_", tag, ".png")), width = 1600, height = 650, res = 110)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), mgp = c(2.2, 0.6, 0))
for (ar in c("visual", "motor")) {
  pm <- prof[area == ar, .(RR = mean(RR), se = sd(RR) / sqrt(.N)), by = .(role, lag_s)]
  pm[is.na(se), se := 0]   # a role with one window
  plot(NA, xlim = range(pm$lag_s), ylim = range(c(pm$RR - 2 * pm$se, pm$RR + 2 * pm$se)), xlab = "lag (s): positive = B follows A",
       ylab = "recurrence rate on the diagonal (%)", main = sprintf("%s: diagonal profile, mean +- 2 SE (%s)", ar,
       main_setting), cex.main = 0.95)
  abline(v = 0, col = "grey80")
  for (ro in intersect(roles, unique(pm$role))) {
    q <- pm[role == ro]
    polygon(c(q$lag_s, rev(q$lag_s)), c(q$RR - 2 * q$se, rev(q$RR + 2 * q$se)),
            col = adjustcolor(role_colour[ro], 0.15), border = NA)
    lines(q$lag_s, q$RR, col = role_colour[ro], lwd = 2)
  }
  legend("topright", intersect(roles, unique(pm$role)), col = role_colour[intersect(roles, unique(pm$role))],
         lwd = 2, bty = "n")
}
dev.off()

# --- 5. Main setting: cross-recurrence plots of one pair --------------------------------------------------
show_s <- 60
f <- folders[grepl(paste0("_", main[1], "$"), folders)][1]
w <- fread(file.path(root, f, "windows.csv"))
env <- fread(file.path(root, f, "envelopes.csv"))
env[, x_a := as.numeric(scale(x_a)), by = window_id][, x_b := as.numeric(scale(x_b)), by = window_id]
ex_pair <- w[, .(n_roles = uniqueN(role)), by = pair][order(-n_roles, pair)]$pair[1]   # a pair with every role
for (ar in c("visual", "motor")) {
  p <- params[folder == f & area == ar & embedding == main[2] & target == as.numeric(main[3])][1]
  ex <- w[pair == ex_pair & area == ar][order(match(role, roles), round, trial), .SD[1], by = role]
  png(file.path(out, sprintf("crp_%s_%s.png", tag, ar)), width = 1400, height = 1450, res = 110)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2.5, 0), mgp = c(2.2, 0.6, 0), pty = "s")
  for (k in seq_len(nrow(ex))) {
    s <- env[window_id == ex$window_id[k]]
    r <- crqa(s$x_a, s$x_b, delay = p$delay, embed = p$embed, radius = p$radius, rescale = 0, normalize = 0,
              mindiagline = 2, minvertline = 2, tw = 0, whiteline = FALSE, recpt = FALSE, side = "both",
              method = "crqa", metric = "euclidean", datatype = "continuous")
    n <- min(nrow(r$RP), round(show_s * p$fs))
    rp <- as.matrix(r$RP[1:n, 1:n])
    t <- (seq_len(n) - 1) / p$fs
    image(t, t, rp * 1, col = c("white", role_colour[ex$role[k]]), useRaster = TRUE,
          xlab = paste("time of", ex$participant_a[k], "(s)"), ylab = paste("time of", ex$participant_b[k], "(s)"),
          main = sprintf("%s: %s r%d t%d\nwhole window: RR %.1f %%, DET %.1f %%", ex$role[k], ex$condition[k],
                         ex$round[k], ex$trial[k], r$RR, r$DET), cex.main = 0.9)
    abline(0, 1, col = "grey40", lty = 3)
  }
  mtext(sprintf("%s, %s %s | %s: delay %d (%.1f s), embedding %d, radius %.2f", ex_pair, ar, ex$band[1], main_setting,
                p$delay, p$delay_s, p$embed, p$radius), outer = TRUE, cex = 0.9)
  dev.off()
}

# --- 6. Outlier windows (main setting) ------------------------------------------------------------------------
# Flagged, not removed: > 3 IQR from the median of their area and role (any measure), or > 20 % marked bad.
d <- results[setting == main_setting & kind == "real"]
flag <- d[, {
  o <- matrix(sapply(measures, function(m) abs(get(m) - median(get(m))) > 3 * IQR(get(m))), nrow = .N)
  .(window_id, outlier_measure = apply(o, 1, function(r) paste(measures[which(r)], collapse = " ")))
}, by = .(session, area, role)]
d <- merge(d, flag, by = c("session", "area", "role", "window_id"))
d[, bad := pmax(bad_fraction_a, bad_fraction_b)]
fwrite(d[outlier_measure != "" | bad > 0.2, .(session, area, role, pair, condition, round, trial, duration_s, RR, DET,
                                             lagRR_2s, bad_fraction_a, bad_fraction_b, outlier_measure)],
       file.path(out, "outliers.csv"))
cat("outlier/bad windows:", d[outlier_measure != "" | bad > 0.2, .N], "of", nrow(d), "\n")
print(summ[contrast == "teammates - opponents" & measure == "lagRR_2s", .(area, setting, n_pairs, mean = signif(mean, 2),
      lo = signif(lo, 2), hi = signif(hi, 2), positive = round(positive, 2))], nrows = 200)
cat("saved to", out, "\n")
