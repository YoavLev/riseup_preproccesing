# Changes: `eeg_data_preprocessing_mne_ordered.ipynb`

Reordered and bug-fixed copy of `eeg_data_preprocessing_mne.ipynb`. The original notebook is unchanged.

New order: 1. load → 2. events → 3. crop + mark breaks → 4. resample → 5. bad channels → 6. ICA → 7. filter + apply ICA → 8. interpolate → 9. re-reference → 10. save + report

## Bug fixes

| # | Where (original) | Change | Why |
|---|---|---|---|
| 1 | ICA, option 1 | Result of `ica.apply` is kept (`raw_clean`) | `raw = ica.apply(raw.copy())` only rebound the loop variable, so the ICA cleaning was lost |
| 2 | ICA, option 1 | `ica.exclude` is one flat list of integers | `[eog_indices, muscle_indices]` was a list of lists |
| 3 | Loading | One participant per run (`SUBJECT`) instead of `range(1, len(files))` | The range skipped the last subject |
| 4 | Loading | No `set_eeg_reference` before `load_data()`; data is preloaded when read | Re-referencing needs preloaded data |
| 5 | Bad channels | `scaling` → `scalings=dict(...)` | `scaling` was undefined (NameError) |
| 6 | Epoching | Tutorial cells removed (`auditory/visual`, `EEG 021`, `pick("mag")`) | Would crash on our data, and ERP epochs don't fit the cRQA analysis |
| 7 | Downsampling | Resample once, on the raw data, in step 4 | The loop resampled the same `epochs` object six times |

## Order and method changes

| # | Where | Change | Why |
|---|---|---|---|
| 8 | Events | Moved right after loading; stored as annotations | Needed for cropping/breaks before ICA; annotations stay correct through cropping and resampling |
| 9 | Events | Each participant's own events, not `raws[0]`'s for everyone | Solo rounds are independent, so the events differ between participants |
| 10 | New: step 2 | Check that all files line up (start time, sampling rate, length, event counts) | Coupling analysis assumes all brains share one time axis |
| 11 | New: step 3 | Crop to first round start → last round end; mark gaps > 15 s between rounds as `BAD_break` | ICA should learn from task data only, not setup or breaks |
| 12 | Downsampling | Moved from the end to step 4 (250 Hz) | Makes ICA and every later step faster |
| 13 | Bad channels | Per participant (`BAD_CHANNELS` dict) instead of `["C4", "TP9"]` for everyone | Bad channels differ between people |
| 14 | Bad channels | Interpolation moved after ICA (step 8) | Interpolated channels add no information and lower the data rank available to ICA |
| 15 | Reference | Mastoids (M1/M2) at load → average reference after interpolation (`REFERENCE` setting) | A noisy mastoid spreads into all channels; ICLabel expects average reference. **To confirm with the lab.** |
| 16 | ICA | One method: ICLabel + EOG correlation as a second opinion + visual check, with a manual override per participant (`ICA_EXCLUDE_MANUAL`) | The original had two alternative options (automatic and corrmap); corrmap refitted ICA and crashed if no EOG component was found |
| 17 | ICA | Picard with `extended=True` (extended Infomax) | ICLabel was trained on extended Infomax; the original markdown said "runica" while the code used plain Picard |
| 18 | ICA | Fitted on a 1–100 Hz, average-referenced copy, EEG channels only, skipping `BAD_break` | What ICLabel expects; EOG/EMG/IMU channels must not enter the decomposition |
| 19 | ICA | Average reference set on the analysis data before applying ICA | ICA must be applied to data with the same reference it was fitted on |
| 20 | Filtering | Notch filter removed | The 40 Hz low-pass already removes 50 Hz line noise |
| 21 | Filtering | Moved before applying ICA; same 0.1–40 Hz | ICA weights come from the separate ICA copy, so the analysis data can be filtered independently |
| 22 | New: step 10 | Saved: cleaned data, ICA object, component table, HTML report, summary CSV (bad channels and excluded components per participant) | Nothing was saved before; these numbers are needed for the Methods |
| 23 | Paths | Outputs go to `<data folder>/derivatives/sub-XX/` | Raw files are never modified, and data stays out of the git repo |
| 24 | All | All parameters in one cell at the top | Every choice in one place, easy to report and change |

## Still to confirm (marked `TODO` in the notebook)
- File names and structure of the `build.py --export-eeg` output: one file per participant, or several people per file?
- Event names and codes from the game repo, in particular the round start/end events.
- EOG channel names (is there an `EOG1`?), and the EMG/IMU channels.
- Reference: average or mastoids?
