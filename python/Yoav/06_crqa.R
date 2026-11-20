# 06_crqa.R: cross-recurrence quantification (cRQA) of every pair x area x trial window of one session,

# Steps:
#   1. z-score each envelope (scale(), base R). optimizeParam() takes its candidate radii from the distances
#      between the values it is given, so they must be on the same scale as what crqa() compares.
#   2. optimizeParam() on every window: delay (auto mutual information), embedding dimension (false nearest
#      neighbours) and the radius that gives 2-5 % recurrence.
#   3. One delay, dimension and radius per area: the median over all windows of that area (all roles pooled,
#      so no role is favoured). Fixed for all windows, so RR can differ between roles.
#   4. crqa() on every window with these parameters; the measures are saved as crqa() returns them.
#
# Usage (on UCloud, from any folder):
#   Rscript 06_crqa.R session03 [cores]
# RT_CRQA_ROOT overrides the folder (e.g. for a test copy).

suppressMessages({
  library(crqa)
  library(data.table)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
session <- if (length(args) >= 1) args[1] else "session03"
cores <- if (length(args) >= 2) as.integer(args[2]) else 8
dir <- file.path(Sys.getenv("RT_CRQA_ROOT", "/work/YoavLevy#3923/derivatives/crqa"), session)
RNGkind("L'Ecuyer-CMRG")  # random numbers that stay reproducible with mclapply
set.seed(97)

# --- Settings -------------------------------------------------------------------------------------------------
# Arguments of crqa(), also used by optimizeParam() (see ?crqa and ?optimizeParam)
normalize <- 0        # the envelopes are z-scored in step 1 already
rescale <- 0          # radius in units of the z-scored envelopes
mindiagline <- 2      # crqa() default
minvertline <- 2      # crqa() default
tw <- 0               # Theiler window: 0 for cross-recurrence (two different people)
# optimizeParam() only
par <- list(method = "crqa", metric = "euclidean", datatype = "continuous",
            maxlag = 50,            # delay: auto mutual information up to 50 samples (2 s at 25 Hz)
            typeami = "mindip",     # ... at its first minimum
            maxEmb = 10,            # embedding dimension: false nearest neighbours for 1..10 dimensions
            fnnpercent = 10,        # ... add dimensions while each removes more than 1/10 of the first share
            radiusspan = 100, radiussample = 40,   # radius search grid (help page: radiusspan e.g. 100)
            normalize = normalize, rescale = rescale, mindiagline = mindiagline, minvertline = minvertline,
            tw = tw, whiteline = FALSE, recpt = FALSE, side = "both")
min_rec <- 2          # radius: recurrence rate between 2 % ...
max_rec <- 5          # ... and 5 % (Coco et al., 2021, examples)

# --- Data -----------------------------------------------------------------------------------------------------
windows <- fread(file.path(dir, "windows.csv"))
envelopes <- fread(file.path(dir, "envelopes.csv"))
envelopes[, x_a := as.numeric(scale(x_a)), by = window_id]   # step 1: z-score per window and person
envelopes[, x_b := as.numeric(scale(x_b)), by = window_id]
series <- split(envelopes, by = "window_id")                   # one table per window, named by window_id
cat(session, ":", nrow(windows), "windows,", cores, "cores\n")

# --- Step 2: optimizeParam() per window -----------------------------------------------------------------------
estimates <- mclapply(windows$window_id, function(id) {
  w <- series[[as.character(id)]]
  p <- optimizeParam(w$x_a, w$x_b, par, min.rec = min_rec, max.rec = max_rec)
  if (is.null(p)) p <- list(radius = NA, emddim = NA, delay = NA)  # no radius found in the 2-5 % range
  data.table(window_id = id, delay = p$delay, embed = p$emddim, radius = p$radius)
}, mc.cores = cores)
estimates <- merge(windows[, .(window_id, area, role)], rbindlist(estimates), by = "window_id")
cat("optimizeParam() found no radius for", sum(is.na(estimates$radius)), "of", nrow(estimates), "windows\n")

# --- Step 3: one set of parameters per area (median over windows) ---------------------------------------------
used <- estimates[, .(delay = round(median(delay, na.rm = TRUE)), embed = round(median(embed, na.rm = TRUE)),
                      radius = median(radius, na.rm = TRUE), n_windows = sum(!is.na(radius))), by = area]
print(used)
fwrite(rbind(estimates[, .(area, window_id, role, delay, embed, radius, used = FALSE)],
             used[, .(area, window_id = NA, role = NA, delay, embed, radius, used = TRUE)]),
       file.path(dir, paste0("crqa_parameters_", session, ".csv")))

# --- Step 4: crqa() per window --------------------------------------------------------------------------------
results <- mclapply(windows$window_id, function(id) {
  w <- series[[as.character(id)]]
  p <- used[area == windows[window_id == id, area]]
  r <- crqa(w$x_a, w$x_b, delay = p$delay, embed = p$embed, radius = p$radius, rescale = rescale,
            normalize = normalize, mindiagline = mindiagline, minvertline = minvertline, tw = tw,
            whiteline = FALSE, recpt = FALSE, side = "both", method = "crqa", metric = "euclidean",
            datatype = "continuous")
  data.table(window_id = id, delay = p$delay, embed = p$embed, radius = p$radius,
             RR = r$RR, DET = r$DET, NRLINE = r$NRLINE, maxL = r$maxL, L = r$L, ENTR = r$ENTR,
             rENTR = r$rENTR, LAM = r$LAM, TT = r$TT)
}, mc.cores = cores)
results <- merge(windows, rbindlist(results), by = "window_id")
fwrite(results, file.path(dir, paste0("crqa_", session, ".csv")))

# --- Overview -------------------------------------------------------------------------------------------------
print(results[, .(n = .N, RR = mean(RR), RR_sd = sd(RR), DET = mean(DET), DET_sd = sd(DET), L = mean(L)),
              by = .(area, role)][order(area, role)], digits = 3)
cat("saved to", dir, "\n")
