# Guide to the preprocessing notebooks

The ICA cleaning follows the **reference pipeline** the lab uses: Stropahl, Bauer, Debener & Bleichner (2018), supplementary material (`notes/papers/eeg_data_sheet.PDF`). Its three MATLAB scripts became three notebooks; two more use ICA as a spatial filter (one occipital alpha component per participant).

| Reference script | Notebook | Run | What it does |
|---|---|---|---|
| `ana01_ICA.m` | `01_ica.ipynb` | per participant | Prepares the data and learns the session ICA |
| `ana02_corrmap.m` | `02_corrmap.ipynb` | once per session | Finds the artifact components with CORRMAP |
| `ana03_preprocessing.m` | `03_clean.ipynb` | per participant | Removes the components and filters the data |
| — | `04a_localizer_ica.ipynb` | per participant | A separate ICA on the visual localizer (instruction video), fitted on the clean data; PNGs for choosing the occipital alpha component |
| — | `04_spatial_filter.ipynb` | per participant | After the choice (`SPATIAL_FILTER_CHOICE`), the alpha envelope of the chosen localizer component at 25 Hz |
| — | `05_windows.ipynb` | once per session | Smooths the envelopes (1 Hz) and downsamples them (5 Hz), then cuts every pair × trial window in session time (with the pair's role) |
| — | `06_crqa.R` + `07_crqa_plots.R` | once, all sessions | R's `crqa` package computes the recurrence measures with one set of settings; then the within-pair contrasts, mixed models and figures |

- **Settings:** every constant is in `config.py`, with its name in the reference scripts in brackets (e.g. `[ana01: LP]`).
- **Helpers (`utils.py`):** `jointprob_keep` (= EEGLAB `pop_jointprob`), `find_bad_channels` (PREP criteria, adapted), `plot_component_maps`.
- **Scripts (`scripts/`):** `setup_job.sh` (Python environment, once per UCloud job); `install_r.sh` (R + crqa + lme4, once per UCloud job); `execute_notebooks.py <notebook> <session> [sub-NNN ...] [--parallel N]` (runs a notebook for everyone in a session); `run_status.py <notebook>` (progress); `suggest_corrmap_templates.py <session>` (ICLabel suggestions for templates); `component_scores.py <session>` + `component_candidates.py <session>` (numbers and shortlist PNGs for choosing `SPATIAL_FILTER_CHOICE`); `final_components.py <session>` (one figure of the chosen components).
- **Open questions** about the adaptation: `ica_reference_questions.md`.
- **Data:** the BIDS dataset on UCloud (read-only), read with `rt_bids`. Participants are `sub-NNN`; their session is `group` in `participants.tsv`.
- **Outputs:** `/work/YoavLevy#3923/derivatives/preproc/sub-NNN/` (`OUT_ROOT`); review PNGs in `derivatives/preproc/spatial_filter_review/<session>/`; executed notebooks and run output in `/work/YoavLevy#3923/logs/`.

```
01_ica (each)                      02_corrmap (session)        03_clean (each)               04a / 04 (each)
load, drop M1/M2 ─► events ─►      templates ─► CORRMAP ─►     apply ICA to 500 Hz data ─►   visual localizer ICA
breaks + crop ─► bad stretches ─►  components to remove        filter 0.1–40 Hz ─►           ─► choose components by eye
bad channels ─► average ref ─►                                 save clean data + report      ─► envelopes
save ─► ICA copy (1–40 Hz, 250 Hz)
─► 1 s epochs, reject > 3 SD ─► Infomax
```

---

## `01_ica.ipynb`: ICA per participant

1. **Load** the BDF with `rt_bids.read_recording` (channels named and typed; the dataset's bad channels in `raw.info["bads"]`), remove M1/M2, set the electrode positions. Result: 27 EEG + 3 EOG channels.
2. **Events:** every event of `events.tsv` is an annotation at its true sample. Counts `acq_gap`s and trials; a check cell compares trial counts across the session.
3. **Breaks and cropping:** gaps > 15 s between trials are labelled `BAD_break` (skipped by later steps), then the data are cropped to the trials (2 s before the first one, the motor baseline). Labels are added *before* cropping, so their times need no correction.
3b. **Bad stretches:** 1 s windows whose median channel SD (1–40 Hz, average reference) is > 50 µV are labelled `BAD_amplitude` (all channels jumping at once, e.g. sub-013). Loss per trial in `<sub>_bad_segments.csv`.
4. **Bad channels:** the dataset's, the automatic ones (`find_bad_channels`) and any in `BAD_CHANNELS` (`config.py`) are **interpolated before ICA** (CORRMAP needs the same channels everywhere). Then the **common average reference**. Prints the median channel SD (units check).
5. **Save** the prepared data (500 Hz) for `03_clean`.
6. **ICA copy:** 1–40 Hz with MNE's `filter()`, then 250 Hz. The high-pass removes slow drift, which would otherwise take many components (question 9).
7. **1 s dummy epochs** (skipping `BAD_` labels); unusual ones (joint probability > 3 SD) are dropped.
8. **Extended Infomax** (= EEGLAB `runica 'extended', 1`). Components = data rank (27 − interpolated − 1). Seed fixed.
9. **Save** the ICA, a PNG of all maps, `_ica_info.csv` (numbers for the Methods); 10. a property plot (map, time course, spectrum) of every component in `components/`.

## `02_corrmap.ipynb`: artifact components, one session at a time

1. Loads the ICA of every participant of the session.
2. **Templates** (`CORRMAP_TEMPLATES` in `config.py`), one clear example per artifact, from that session; suggestions from `suggest_corrmap_templates.py`:
   - **eye blink:** strong over the forehead, both sides the same sign;
   - **lateral eye movement:** front left vs. front right, opposite signs;
   - **heartbeat:** left out so far (no clear component); **no muscle template** (supervisor: it can remove brain signal).
   MNE counts components from **0** (EEGLAB from 1).
3. **CORRMAP** finds components with a similar map (correlation ≥ 0.8) in every participant, in two rounds.
4. **Save** `_badcomps.csv` per participant (matches + `EXTRA_EXCLUDE` − `EXTRA_KEEP`) and the maps PNG with the removed components in red.

## `03_clean.ipynb`: cleaning per participant

1. Loads the prepared data, the ICA and the components to remove.
2. **Removes the components** from the 500 Hz data (the ICA weights are a fixed spatial unmixing, so they apply to the original data). Overlay plot: blinks should disappear, the rest barely change.
3. **Filters** 0.1–40 Hz with MNE's `filter()`.
4. **Saves** `_clean_raw.fif` (input for everything after), an HTML report, and a row in `preprocessing_summary.csv`.

## `04a_localizer_ica.ipynb`: ICA on the visual localizer

Only **occipital alpha** is analysed (no motor components). Per participant, on the clean data, one ICA on the **visual localizer** (`LOCALIZER` in `config.py`): the instruction video, from its first frame on the screen (photodiode; session04 has none, so the phase's start marker). Same settings as `01` (1–40 Hz, 250 Hz, 1 s epochs, 3 SD, extended Infomax); `BAD_amplitude` skipped, `BAD_break` ignored (the video lies in a break). Components = rank of the clean data (`01` components − removed). Saves `<sub>_visual_localizer_ica.fif`, `<sub>_localizer_ica_info.csv` and PNGs in `spatial_filter_review/<session>/localizer/`. The filter is learned from data the analysis does not test.

**Choosing the component:** `scripts/component_scores.py <session>` (alpha peak above 1/f, 20–40 Hz slope, posterior share of the map, ICLabel) and `scripts/component_candidates.py <session>` (shortlist: strongest electrode O1/O2, ≥ 35 % of the map posterior, brain-like slope, alpha peak ≥ 1 dB; constants in `config.py`), then by eye. A participant without a clear occipital alpha component is excluded.

## `04_spatial_filter.ipynb`: the occipital alpha envelope

After the choice is written into `SPATIAL_FILTER_CHOICE` (`{participant: component of the visual localizer ICA}`), takes the component's time course on the clean data, band-passes it (`ALPHA_BAND`, 8–12 Hz, 1 Hz slopes), takes its Hilbert envelope on the continuous data and downsamples it to `ENVELOPE_SFREQ` (25 Hz, so cRQA stays feasible). Saved as `<sub>_visual_alpha_envelope_raw.fif`, with all events. Sections 2b and 2c only explain the steps (and why the envelope is smoothed in `05`).

## `05_windows.ipynb`, `06_crqa.R`, `07_crqa_plots.R`: cross-recurrence per pair

`05_windows` (Python, because the windows are cut with `rt_bids`):
1. Loads the alpha envelope of every participant with a choice (checks that `04` used what `config.py` says); pairs = any two of them.
2. Trials from `rt_bids.session_events` (session time); per pair and trial the window both were playing; role: `teammates` / `opponents` (same / different `team_id` in a joint trial), `solo`, `baseline`.
3. Session time → envelope sample with `rt_bids.from_session_time` (drift and gaps corrected) minus the clean data's `first_samp`; not `crop_session` (it ignores `first_samp`). The share of bad time per window is recorded.
4. Smooths each envelope (`CRQA_SMOOTH_HZ` = 1 Hz, zero-phase FIR) and downsamples it (`CRQA_SFREQ` = 5 Hz); saves `derivatives/crqa/<session>_smooth1hz_5hz/windows.csv` and `envelopes.csv`.

`06_crqa.R` (R, `crqa` package; once per job `bash scripts/install_r.sh`, then `Rscript 06_crqa.R session02,session03,... [cores]`), all sessions together, one set of settings:
1. z-scores each envelope per window; embedding m = 1 (no embedding, so the delay plays no role).
2. One radius (per area; only `visual` now): the median over all windows of the 5 % quantile of the A–B distances (target ~5 % recurrence), fixed for all roles.
3. Minimum line 0.4 s. Measures: `lagRR_2s` (recurrence within ±2 s of the same moment, primary), `RR` of the whole plot (control), DET, L, LAM; surrogate pairs (the same pair in another block, same trial number).
4. Saves `derivatives/crqa/final_smooth1hz_5hz/` (`parameters.csv`, `windows.csv`, `results.csv`, `profiles.csv`).

`07_crqa_plots.R`: within-pair contrasts (teammates − opponents, teammates − solo, joint − surrogate), mixed models (lme4) and figures, in the same folder.

## Order of work

1. `01_ica` for every participant of a session (check bad channels and bad stretches).
2. Choose the templates, run `02_corrmap`; review the red components, add misses to `EXTRA_EXCLUDE` and false matches to `EXTRA_KEEP`, rerun.
3. `03_clean` for every participant.
4. `04a_localizer_ica` for everyone; `component_scores.py` + `component_candidates.py`; choose by eye (`SPATIAL_FILTER_CHOICE`); `04_spatial_filter` (envelopes).
5. `05_windows` once per session; then `06_crqa.R` over all sessions and `07_crqa_plots.R`.
