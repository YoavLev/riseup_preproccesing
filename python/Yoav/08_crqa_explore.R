# 08_crqa_explore.R: cRQA of one window folder of 05_windows (smoothed and downsampled envelopes) with parameters
# chosen ONCE per area from all windows pooled, for two target recurrence rates, plus surrogate pairs.
# Written for the parameter exploration of 2026-11-24 (compare smoothing levels; 06_crqa.R stays as it was).
#
# Why pooled instead of 06_crqa.R's per-window optimizeParam() + median: the per-window estimates split into two
# groups (embedding 1 with radius ~0.05, embedding 5-6 with radius ~1.3-1.6), and the median of each parameter
# taken separately mixed the groups. Here each parameter is taken from one average curve, in order, each step using
# the previous one (the order of Coco et al., 2021: delay -> embedding -> radius).
#
# Steps (per area, all roles and both people pooled, so no role is favoured):
#   1. z-score each envelope per window (as 06_crqa.R).
#   2. Delay: average mutual information (AMI) of each envelope with itself, averaged over all windows and both
#      people; delay = the FIRST local minimum of the average curve (Coco et al., 2021; "first minimum" rule of
#      Fraser & Swinney, 1986, unverified - check before citing). At this lag the envelope's value says least about
#      its earlier value, so the embedded coordinates add new information.
#   3. Embedding dimension, two versions (both are computed and compared):
#      Both from false nearest neighbours (tseriesChaos::false.nearest with the settings optimizeParam() uses:
#      rt = 10, eps = SD / 10, t = 0) with that delay, for 1..max_emb dimensions, averaged the same way, and the rule
#      of optimizeParam() (fnnpercent = 10; Coco et al., 2021): the dimension after the last drop of the curve larger
#      than 1/10 of its value at dimension 1. They differ only in which dimensions count:
#      "fnn":       only dimensions where most windows still have neighbours (in higher dimensions eps = SD/10 leaves
#                   no neighbours, false.nearest() gives NA, and the average rests on a few windows);
#      "crqa_rule": every dimension with a value in any window (what optimizeParam() does with one window).
#      In a test (2026-11-24) the FNN share stayed at 80-90 % instead of falling towards 0: the envelopes behave
#      like noise at this time scale (FNN is meant for nearly deterministic signals), so "fnn" gave 1 (no embedding,
#      the usual choice for noisy data; Iwanski & Bradley, 1998, unverified - check before citing) and "crqa_rule"
#      a larger dimension from the last, thinly-supported drop.
#   4. Radius (for each embedding version): for each window, the distance between A's and B's embedded states below which target % of all
#      state pairs fall (the target-quantile of the cross-distance matrix); radius = median over the windows of the
#      area. Fixed for all windows, so %REC differs between windows and roles (it is a dependent variable).
#      Two targets: 2 % and 5 % (the range of Coco et al., 2021).
#   5. crqa() of every window with these parameters: RR (%REC), DET, L, ... ; minimum line length 2 samples
#      (crqa default; at a sampling rate of 5 x the smoothing cut-off, 2 samples = 0.4 / cut-off s, the time scale
#      the smoothing leaves). From the recurrence matrix also: the diagonal profile (recurrence rate per lag, -10..10
#      s) and %REC within +-2 s of the same moment (lagRR_2s: coupling close in time, not any time in the trial).
#   6. Surrogate pairs (joint trials only): A's envelope in one joint trial with B's envelope in a joint trial of
#      another block (same trial number), cut to the shorter. Same people, same kind of task, not the same moment:
#      what each person's own dynamics produce. Real joint pairs minus surrogates = coupling of the same moment.
#
# Usage (on UCloud, after scripts/install_r.sh):  Rscript 08_crqa_explore.R session03_smooth1hz_5hz [cores]
# Output: derivatives/crqa/<folder>/explore/  parameters.csv, curves.csv, results.csv, profiles.csv

suppressMessages({
  library(crqa)
  library(tseriesChaos)
  library(data.table)
  library(parallel)
  library(Matrix)
})
setDTthreads(1)   # data.table's OpenMP threads make forked mclapply() workers hang (see 06_crqa.R)

args <- commandArgs(trailingOnly = TRUE)
folder <- if (length(args) >= 1) args[1] else "session03_smooth1hz_5hz"
cores <- if (length(args) >= 2) as.integer(args[2]) else 8
dir <- file.path(Sys.getenv("RT_CRQA_ROOT", "/work/YoavLevy#3923/derivatives/crqa"), folder)
out <- file.path(dir, "explore")
dir.create(out, showWarnings = FALSE)
RNGkind("L'Ecuyer-CMRG")
set.seed(97)

# --- Settings -------------------------------------------------------------------------------------------------
max_lag_s <- 5        # AMI up to 5 s
max_emb <- 8          # FNN for 1..8 dimensions
fnnpercent <- 10      # embedding: last drop larger than 1/10 of the FNN share at dimension 1 (optimizeParam default)
targets <- c(2, 5)    # target recurrence rates [%] for the radius
mindiagline <- 2      # minimum diagonal line [samples] (crqa default)
minvertline <- 2
band_s <- 2           # lagRR: recurrence within +-2 s of the same moment
profile_s <- 10       # diagonal profile saved for -10..10 s

# --- Data -----------------------------------------------------------------------------------------------------
windows <- fread(file.path(dir, "windows.csv"))
envelopes <- fread(file.path(dir, "envelopes.csv"))
envelopes[, x_a := as.numeric(scale(x_a)), by = window_id]   # step 1
envelopes[, x_b := as.numeric(scale(x_b)), by = window_id]
series <- split(envelopes, by = "window_id")
fs <- round(median(windows$n_samples / windows$duration_s), 3)
cat(folder, ":", nrow(windows), "windows at", fs, "Hz,", cores, "cores\n")
t0 <- Sys.time()
elapsed <- function() sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

# --- Step 2: delay from the average AMI curve -----------------------------------------------------------------
max_lag <- round(max_lag_s * fs)
ami <- mclapply(windows$window_id, function(id) {
  w <- series[[as.character(id)]]
  rbind(as.numeric(mutual(w$x_a, lag.max = max_lag, plot = FALSE)),
        as.numeric(mutual(w$x_b, lag.max = max_lag, plot = FALSE)))
}, mc.cores = cores)
first_min <- function(y) {   # first lag (from 1) where the curve stops falling; NA if it falls all the way
  k <- which(diff(y) >= 0)[1]
  if (is.na(k)) NA_integer_ else as.integer(k - 1)   # y[1] is lag 0
}
params <- list(); curves <- list()
for (ar in unique(windows$area)) {
  ids <- which(windows$area == ar)
  m <- colMeans(do.call(rbind, ami[ids]))
  d <- first_min(m)
  if (is.na(d)) d <- which(m < m[1] / exp(1))[1] - 1   # fallback: first lag below 1/e of lag 0 (reported)
  params[[ar]] <- list(area = ar, delay = d, delay_s = d / fs)
  curves[[length(curves) + 1]] <- data.table(area = ar, curve = "AMI", x = 0:max_lag, x_s = (0:max_lag) / fs, y = m)
}
cat("delays:", paste(sapply(params, function(p) sprintf("%s %d (%.2f s)", p$area, p$delay, p$delay_s)), collapse = ", "),
    "|", elapsed(), "\n")

# --- Step 3: embedding from the average FNN curve -------------------------------------------------------------
fnn <- mclapply(windows$window_id, function(id) {
  w <- series[[as.character(id)]]
  d <- params[[windows[window_id == id, area]]]$delay
  f <- function(x) 100 * as.numeric(false.nearest(x, m = max_emb, d = d, t = 0, rt = 10, eps = sd(x) / 10)[1, ])
  rbind(f(w$x_a), f(w$x_b))
}, mc.cores = cores)
for (ar in names(params)) {
  ids <- which(windows$area == ar)
  all <- do.call(rbind, fnn[ids])
  m <- colMeans(all, na.rm = TRUE)
  rule <- function(y) {                                  # optimizeParam(): dimension after the last large drop
    y <- y[!is.na(y)]
    drops <- which(diff(y) < -y[1] / fnnpercent)         # drop from dimension k to k + 1
    if (length(drops)) max(drops) + 1L else 1L
  }
  params[[ar]]$embed_crqa_rule <- rule(m)
  params[[ar]]$embed_fnn <- rule(ifelse(colMeans(is.na(all)) > 0.5, NA, m))   # most windows without neighbours: out
  curves[[length(curves) + 1]] <- data.table(area = ar, curve = "FNN_share_NA", x = 1:max_emb, x_s = NA,
                                             y = 100 * colMeans(is.na(all)))
  curves[[length(curves) + 1]] <- data.table(area = ar, curve = "FNN", x = 1:max_emb, x_s = NA, y = m)
}
cat("embedding (fnn / crqa_rule):", paste(sapply(params, function(p) sprintf("%s %d / %d", p$area, p$embed_fnn,
    p$embed_crqa_rule)), collapse = ", "), "|",
    elapsed(), "\n")

# --- Step 4: radius for each target recurrence rate -----------------------------------------------------------
embed_ts <- function(x, m, d) {   # rows = states (x[t], x[t+d], ...), as crqa(); embedd() fails for m = 1
  if (m == 1) matrix(x, ncol = 1) else embedd(x, m = m, d = d)
}
cross_dist <- function(a, b) {                           # Euclidean distance of every state of A to every state of B
  d2 <- outer(rowSums(a^2), rowSums(b^2), "+") - 2 * tcrossprod(a, b)
  sqrt(pmax(d2, 0))
}
embeddings <- c("fnn", "crqa_rule")
embed_of <- function(p, version) as.integer(p[[paste0("embed_", version)]])
radii <- mclapply(windows$window_id, function(id) {
  w <- series[[as.character(id)]]
  p <- params[[windows[window_id == id, area]]]
  unlist(lapply(embeddings, function(v) {
    m <- embed_of(p, v)
    D <- cross_dist(embed_ts(w$x_a, m, p$delay), embed_ts(w$x_b, m, p$delay))
    quantile(D, targets / 100, names = FALSE)
  }))
}, mc.cores = cores)
radii <- do.call(rbind, radii)   # columns: fnn 2 %, fnn 5 %, crqa_rule 2 %, crqa_rule 5 %
settings <- CJ(area = names(params), embedding = embeddings, target = targets, sorted = FALSE)
settings[, `:=`(delay = sapply(area, function(a) params[[a]]$delay),
                delay_s = sapply(area, function(a) params[[a]]$delay_s),
                embed = mapply(function(a, v) embed_of(params[[a]], v), area, embedding),
                radius = mapply(function(a, v, tg) median(radii[windows$area == a,
                                (match(v, embeddings) - 1) * length(targets) + match(tg, targets)]),
                                area, embedding, target))]
parameters <- settings
parameters[, `:=`(folder = folder, fs = fs, mindiagline = mindiagline, mindiagline_s = mindiagline / fs)]
print(parameters)
fwrite(parameters, file.path(out, "parameters.csv"))
fwrite(rbindlist(curves), file.path(out, "curves.csv"))
cat("radii done |", elapsed(), "\n")

# --- Steps 5 and 6: crqa() of every window and of the surrogate pairs -----------------------------------------
band <- round(band_s * fs)
max_prof <- round(profile_s * fs)
measure <- function(xa, xb, p) {
  r <- crqa(xa, xb, delay = p$delay, embed = p$embed, radius = p$radius, rescale = 0, normalize = 0,
            mindiagline = mindiagline, minvertline = minvertline, tw = 0, whiteline = FALSE, recpt = FALSE,
            side = "both", method = "crqa", metric = "euclidean", datatype = "continuous")
  n <- nrow(r$RP)
  nz <- summary(r$RP)                       # recurrent points: i = A's time, j = B's time
  lag <- nz$j - nz$i                        # positive = B's state matches A's earlier state (B follows A)
  counts <- tabulate(lag[abs(lag) <= max_prof] + max_prof + 1, nbins = 2 * max_prof + 1)
  cells <- n - abs(-max_prof:max_prof)      # cells on each diagonal
  in_band <- abs(-max_prof:max_prof) <= band
  list(values = data.table(RR = r$RR, DET = r$DET, L = r$L, L_s = r$L / fs, maxL = r$maxL, ENTR = r$ENTR,
                           LAM = r$LAM, TT = r$TT, lagRR_2s = 100 * sum(counts[in_band]) / sum(cells[in_band])),
       profile = 100 * counts / cells)
}

# Surrogate partners: for joint trial (round r, trial t) of a pair, B's envelope from (other round, same t)
joint <- windows[condition == "joint", .(window_id, area, pair, round, trial)]
surrogates <- merge(joint, joint, by = c("area", "pair", "trial"), suffixes = c("_a", "_b"), allow.cartesian = TRUE)
surrogates <- surrogates[round_a != round_b]
jobs <- rbind(windows[, .(kind = "real", window_id, window_b = window_id, area)],
              surrogates[, .(kind = "surrogate", window_id = window_id_a, window_b = window_id_b, area)])
cat(nrow(jobs), "cRQA jobs (", sum(jobs$kind == "real"), "real,", sum(jobs$kind == "surrogate"), "surrogate ) x",
    length(embeddings) * length(targets), "settings\n")

computed <- mclapply(seq_len(nrow(jobs)), function(k) {
  j <- jobs[k]
  xa <- series[[as.character(j$window_id)]]$x_a
  xb <- series[[as.character(j$window_b)]]$x_b
  n <- min(length(xa), length(xb))
  lapply(which(settings$area == j$area), function(s) {
    p <- settings[s]
    m <- measure(xa[1:n], xb[1:n], p)
    list(values = cbind(j, embedding = p$embedding, target = p$target, n_samples = n, m$values),
         profile = if (j$kind == "real") data.table(window_id = j$window_id, embedding = p$embedding,
                                                     target = p$target, lag_s = (-max_prof:max_prof) / fs,
                                                     RR = m$profile))
  })
}, mc.cores = cores)
failed <- sapply(computed, inherits, "try-error")
if (any(failed)) stop(sum(failed), " cRQA jobs failed, e.g. ", computed[[which(failed)[1]]])
computed <- unlist(computed, recursive = FALSE)
results <- rbindlist(lapply(computed, `[[`, "values"))
results <- merge(results, windows[, .(window_id, session, pair, participant_a, participant_b, condition, round,
                                      trial, role, duration_s, bad_fraction_a, bad_fraction_b)], by = "window_id")
results[kind == "surrogate", role := "surrogate"]
fwrite(results, file.path(out, "results.csv"))
fwrite(rbindlist(lapply(computed, `[[`, "profile")), file.path(out, "profiles.csv"))

print(results[, .(n = .N, RR = median(RR), DET = median(DET), L_s = median(L_s), lagRR_2s = median(lagRR_2s)),
              by = .(area, embedding, target, role)][order(area, embedding, target, role)], digits = 3)
cat("saved to", out, "|", elapsed(), "\n")
