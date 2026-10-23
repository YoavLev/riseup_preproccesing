# Guide to `eeg_data_preprocessing_mne_ordered.ipynb`

What each cell does and why. Headings follow the notebook's own headings; "Cell N" counts all cells from 0 (markdown included), in the order they appear.

**The big picture.** The notebook takes the raw EEG of **one participant** and produces a **clean, continuous recording** of the task, plus a record of every cleaning decision. It does not compute any coupling yet: that happens in a separate analysis notebook. Run it once per participant.

```
raw .fif ─► 1 load ─► 2 events ─► 3 crop + breaks ─► 4 resample ─► 5 bad channels
        ─► 6 ICA (on a copy) ─► 7 filter + apply ICA ─► 8 interpolate ─► 9 re-reference ─► 10 save
```

---

## Introduction (cells 0–4)

**Cells 0, 1, 3 (markdown):** title, input files, how to use the notebook, the flow.

**Cell 2: build.py instructions.** Not real code: a text block (a string) reminding you how the input files are made in the RiseTogether-analysis repo. `build.py --export-eeg` writes three files per participant; this notebook uses only the `_eeg_raw.fif` one.

**Cell 4: imports.**
- `mne`: the EEG library. `numpy` (`np`): arrays and math. `pandas` (`pd`): tables.
- `ICA` and `label_components` (ICLabel): for step 6.
- `%matplotlib inline`: plots appear inside the notebook. Change it to `%matplotlib qt` to get interactive plot windows (e.g., to click bad channels).

---

## 0. Parameters (cell 6)

**Purpose:** every setting in one place. Nothing else in the notebook should contain a "magic number". When you write the Methods, this cell is your checklist.

| Group | What it sets |
|---|---|
| Participant | `SESSION` and `IPAD` pick the file; `SUBJECT` (e.g. `s1_pia`) names all outputs |
| Paths | Where the input files are (`out/eeg/`), the channel map (`cap_layout.csv`), and where outputs go (`out/derivatives/sub-<SUBJECT>/`). `mkdir` creates the output folder if needed |
| Channels | The electrode-position template (`MONTAGE`), which channels are EOG |
| `BAD_CHANNELS` | **Filled in by hand** after step 5, one line per participant |
| Events | Which conditions to keep, the start/end event names, which annotations are sync-only, the minimum break length |
| Filtering | Resampling rate (250 Hz); the ICA copy's filter (1–100 Hz); the analysis filter (0.1–40 Hz); the final reference |
| ICA | How many components (95% of the variance), the artifact threshold for fitting, which ICLabel classes count as artifacts and how sure ICLabel must be (0.8), a manual override per participant, the random seed |

**Concept: why a random seed?** ICA starts from a random guess and improves it step by step. A fixed seed (`RANDOM_STATE = 97`) makes the starting guess the same every time, so rerunning the notebook gives exactly the same components.

---

## 1. Loading the raw data (cell 8)

**Purpose:** read this participant's file and make the channels usable.

1. **Find the file.** `glob` searches `out/eeg/` for `s<SESSION>_<IPAD>_*_eeg_raw.fif` (the `*` is the recording name, not known in advance). `assert` stops the notebook with a clear message if there isn't exactly one match.
2. **Load it.** `preload=True` reads all the data into memory, which filtering and ICA need.
3. **Rename the channels.** The file uses hardware names (`p0_ch1`, `p0_ch2`, …). `cap_layout.csv` maps them to electrode names (`O1`, `Pz`, …), and `rename_channels` applies that map.
4. **Set the channel types.** The file marks all 32 channels as EEG. EOG2–4 are eye electrodes, so they are re-typed as `eog`; this keeps them out of ICA and the average reference.
5. **Set the montage.** Gives every electrode a 3D position on the head. Needed for topography plots, ICLabel and interpolation.

**Check:** the printout shows 32 channels at 500 Hz; no error about unknown channel names.

---

## 2. Events (cells 9–12)

**Concept: annotations.** An annotation is a label attached to a moment in the recording (onset, duration, text). `build.py` already placed every game event at the correct sample of this file, for example `TRIAL_START/joint/r1/t2`. Because they are stored with the data, they stay correct when the data is cropped or resampled.

**Cell 10: check and clean the events.**
1. `assert`: a file with no annotations comes from a headset `build.py` could not align. Its data can't be used, so the notebook stops.
2. Remove the **sync annotations** (`trig_in*`, `GPIO_TRIGGER/`, `PHOTON/`). `build.py` used them to align the EEG with the game; they aren't events for the analysis. `str.match` finds them with a regular expression; `delete` removes them.
3. Count the trials per condition (the text between the first and second `/` of each `TRIAL_START` label).
4. Plot the trial starts and ends over time.

**Check:** the trial counts match the QA report; the plot shows the conditions in the expected order (motor baseline → individual → joint).

**Cell 12: completeness check across participants.** Opens every file of the session *without loading the data* (`preload=False`, fast) and builds a table: file, sampling rate, number of channels, aligned or not, trials per condition. All six rows should agree. It deliberately does **not** compare start times: the headsets were switched on one at a time, so their files start at different moments.

---

## 3. Cropping to the task and marking breaks (cells 13–14)

**Purpose:** ICA should learn only from task data. The recording covers the whole session (welcome, surveys, instruction video, debriefing), where people talk and move in ways that don't happen during play.

**Cell 14, step by step:**
1. **`trial_times`** (a small helper function): returns the times, in seconds from the start of the data, of all `TRIAL_START/` (or `TRIAL_END/`) events of the conditions in `CONDITIONS`.
2. **Find the breaks.** For each trial end, look at the next trial start. If the gap is longer than `BREAK_MIN_DURATION` (15 s), it's a break: a survey, results screen, instructions or team reassignment. A 1 s margin is kept on each side.
3. **Mark them as `BAD_break`.** Annotations starting with `BAD` are skipped automatically by ICA and many later MNE functions. The data itself isn't deleted: the label just says "ignore this part".
4. **Crop** to the task: from 2 s before the first trial start to 2 s after the last trial end.

**Why breaks are marked before cropping:** the files have no recording date (`meas_date`). In that case, MNE shifts all events by the cropped amount if annotations are added after a crop. Marking first and cropping second avoids that bug.

**Check:** the printed break total is plausible; in later plots, `BAD_break` regions sit between trials, never inside one.

---

## 4. Downsampling (cell 16)

`raw.resample(250)` keeps 250 samples per second instead of 500. Alpha and beta (up to 30 Hz) need far less: a signal can be represented up to half the sampling rate (125 Hz here). Half the samples makes every later step about twice as fast. MNE low-pass filters first to avoid aliasing (fast signals masquerading as slow ones), and the annotations are kept.

---

## 5. Identifying bad channels (cells 17–19)

**Concept: a bad channel** is an electrode that recorded mostly noise for much of the session (poor contact, broken wire). Unlike a blink, which shows up in all frontal channels together, a bad channel misbehaves **on its own**.

**Cell 18: inspect.**
- `raw.plot`: scroll through the signals. Look for flat lines, constant large noise, or drifts and jumps that only one channel has.
- `compute_psd().plot()`: the power spectrum per channel. A bad channel is a line far above or below all others, or with a sharp 50 Hz peak the others don't have.

**Cell 19: apply.** Copies this participant's list from `BAD_CHANNELS` into `raw.info["bads"]`. From then on MNE knows these channels are bad: ICA skips them and step 8 replaces them.

**Workflow:** run cell 18 → write the channels into `BAD_CHANNELS` in the parameters cell → re-run the parameters cell → run cell 19. (If you clicked channels in an interactive plot, copy them into the dictionary, because cell 19 overwrites the clicks.)

---

## 6. ICA (cells 20–24)

**Concept: ICA (independent component analysis).** Each electrode records a *mixture* of sources: brain activity, eye movements, muscles, heartbeat, line noise. ICA unmixes the channels into **independent components**, each with a scalp map (where it is strongest) and a time course (what it does). Artifact components can then be removed and the remaining components mixed back into channels.

**Cell 21: fit ICA, on a copy.**
1. **Make a copy** filtered at 1–100 Hz and set to average reference. This copy is only for learning the unmixing:
   - The **1 Hz high-pass** removes slow drifts, which carry a lot of power and would otherwise use up components.
   - Keeping frequencies **up to 100 Hz** helps identify muscle components (they show broad high-frequency power).
   - 1–100 Hz with **average reference** is what ICLabel was trained on.
2. **Set up ICA:**
   - `n_components=0.95`: keep as many components as needed to explain 95% of the variance.
   - `method="infomax", extended=True`: extended Infomax, the algorithm ICLabel was trained on. "Extended" lets it separate both peaky sources (blinks) and flatter ones (line noise).
   - `random_state`: the fixed seed (see Parameters).
3. **Fit:** `picks="eeg"` uses EEG channels only (not EOG, not bad channels); `reject_by_annotation=True` skips the `BAD_break` parts; `reject` skips 2 s pieces with a peak-to-peak range above 200 µV (gross artifacts like cable pulls).

**Cell 22: label the components automatically.**
1. **ICLabel** gives each component a label (brain, eye blink, muscle artifact, heart beat, line noise, channel noise, other) and a probability.
2. A component is flagged if its label is an artifact class **and** the probability is at least 0.8.
3. **Second opinion:** `find_bads_eog` flags components whose time course correlates strongly with the EOG channels (eye activity).
4. The union of both lists is `auto_exclude`. The table `components` shows everything.

**Cell 23: look at them.** Topographies of all components, their time courses, and detailed plots of the flagged ones. What to expect:
- **Blink:** strong frontal scalp map; sharp spikes in time.
- **Horizontal eye movement:** left–right opposite poles at the front; step-like time course.
- **Muscle:** map at the edge of the head; high, flat-to-rising power at high frequencies.
- **Brain:** smooth dipolar map; a clear alpha peak around 10 Hz in the spectrum. Keep these.

**Cell 24: decide and save.**
- If `ICA_EXCLUDE_MANUAL` has a list for this participant, that list wins; otherwise `auto_exclude` is used. This is where your visual judgment overrides the automation.
- Saves the ICA solution (`_ica.fif`) and the component table with an `excluded` column (`_ica_components.csv`) for the Methods.

---

## 7. Filtering the analysis data and applying ICA (cells 25–26)

**Purpose:** clean the data you will actually analyze.
1. **Filter a new copy at 0.1–40 Hz.** This is the analysis data, separate from the ICA copy. The 40 Hz low-pass also removes the 50 Hz line noise, so no notch filter is needed.
2. **Set the average reference**, because ICA was fitted on average-referenced data. The unmixing only works on data with the same reference.
3. **Apply ICA:** removes the excluded components and mixes the rest back into channels.
4. **Overlay plot:** the signal before and after cleaning. Blinks should disappear; everything else should look almost unchanged.

**Concept: why fit on one copy and apply to another?** The ICA copy is filtered for *learning* the sources well; the analysis copy is filtered for *your analysis*. ICA is a fixed linear unmixing, so weights learned on one can be applied to the other.

---

## 8. Interpolating bad channels (cell 28)

Saves the list of bad channels (for the summary), then **interpolates** them: each bad channel is replaced by an estimate from its neighbors, weighted by distance (this is why the montage was needed). Afterwards every participant has the same 32-channel layout, which simplifies ROIs and comparisons. It's done after ICA because interpolated channels contain no new information and would only make ICA's job harder.

---

## 9. Re-referencing (cell 30)

**Concept: the reference.** EEG measures voltage *differences*, so every channel is relative to something. With the **average reference**, each channel is expressed relative to the average of all EEG channels. Done after interpolation so the repaired channels are part of the average. The alternative (`["M1", "M2"]`, the mastoids) is still a TODO to confirm with the lab.

The plot is a final visual check of the clean data.

---

## 10. Saving the cleaned data and a quality report (cells 31–32)

Writes three things to `out/derivatives/sub-<SUBJECT>/`:
1. **`_clean_raw.fif`:** the cleaned continuous data, with all events still attached. This is the input for the analysis notebook.
2. **`_report.html`:** an MNE report with the data before and after cleaning (including spectra) and the ICA components. Open it in a browser; it can serve as a supplementary figure.
3. **`preprocessing_summary.csv`** (one folder up, shared by all participants): one row per participant with bad channels, number of components, excluded components, break time and duration. Re-running a participant replaces its row. These are the numbers you report in the Methods (e.g. "on average 1.2 channels were interpolated and 6 components removed").

**Check:** in the report, the spectrum after cleaning still has an alpha peak around 10 Hz. If it disappeared, ICA removed brain signal: look again at the excluded components.

---

## Next: analysis (cell 33)

A reminder of what the separate analysis notebook will do: alpha and beta envelopes on the continuous data, cut into trials, aligned between participants on the same events (never on file start), teammate/opponent labels per trial from the `trials` table, then cRQA.
