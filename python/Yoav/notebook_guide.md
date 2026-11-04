# Guide to the preprocessing notebooks

The preprocessing follows the **reference pipeline** the lab uses for ICA: Stropahl, Bauer, Debener & Bleichner (2018), supplementary material (`notes/papers/eeg_data_sheet.PDF`). Its three MATLAB scripts became three notebooks:

| Reference script | Notebook | Run | What it does |
|---|---|---|---|
| `ana01_ICA.m` | `01_ica.ipynb` | once **per participant** | Prepares the data and learns the ICA |
| `ana02_corrmap.m` | `02_corrmap.ipynb` | once, for **all participants** together | Finds the artifact components with CORRMAP |
| `ana03_preprocessing.m` | `03_clean.ipynb` | once **per participant** | Removes the components and filters the data |
| — (next step) | `04_spatial_filter.ipynb` | once **per participant** | Uses ICA as a spatial filter: saves the candidate visual (alpha) and motor (mu/beta) components as PNG for the review (`spatial_filter_review/<session>/`); after the choice (`SPATIAL_FILTER_CHOICE`), computes their envelopes |

- **Settings:** every constant is in `config.py`, with the name it has in the reference scripts in brackets (e.g. `[ana01: LP]`).
- **Helpers:** `utils.py` has Python versions of two EEGLAB functions: `firws_filter` (= `pop_firws`) and `jointprob_keep` (= `pop_jointprob`).
- **Open questions** about the adaptation: `ica_reference_questions.md`.
- **Data:** the BIDS dataset on UCloud (`/work/RiseTogetherEEG/VideoGameParadigm-EEG` or `/work/VideoGameParadigm-EEG`, depending on how it is mounted; read-only), read with `rt_bids` (the dataset's `code/rt_bids.py`). Participants are `sub-NNN`; their session is `group` in `participants.tsv`.
- **Outputs:** `/work/YoavLevy#3923/derivatives/preproc/sub-NNN/` (`OUT_ROOT` in `config.py`); executed notebooks and script output in `/work/YoavLevy#3923/logs/`.
- **Scripts** (`scripts/`): `setup_job.sh` (environment, once per UCloud job), `execute_notebooks.py <notebook> <session> [sub-NNN ...]` (runs a per-participant notebook for everyone in a session).
- The earlier single notebook `eeg_data_preprocessing_mne_ordered.ipynb` is **superseded** by these three.

```
01_ica (each participant)           02_corrmap (all)            03_clean (each participant)
load, drop M1/M2 ──► events ──►     choose templates ──►        apply ICA to 500 Hz data ──►
breaks + crop ──► bad channels ──►  CORRMAP (r ≥ 0.8) ──►       low-pass 40 / high-pass 0.1 ──►
save prepared data ──► ICA copy     components to remove        save clean data + report
(1–40 Hz, 250 Hz) ──► 1 s epochs,
reject > 3 SD ──► Infomax ICA
```

---

## `01_ica.ipynb`: ICA per participant

**0. Participant.** `SUBJECT` (e.g. `sub-001`) names all outputs. Set it in the cell, or through the environment variable `RT_SUBJECT` (what `execute_notebooks.py` does).

**1. Loading.** Reads the participant's BDF with `rt_bids.read_recording` (channels already named and typed; channels marked bad in `channels.tsv` already in `raw.info["bads"]`), loads it into memory, **removes the mastoids M1 and M2**, makes sure EOG2–4 are eye channels (so they stay out of the ICA), and sets the electrode positions. Result: 27 EEG + 3 EOG channels. It also tags the data with the path of the recording, so later steps can convert sample positions to session time.

**2. Events.** Every event of the recording's `events.tsv` is an *annotation* (a label at a moment in the recording), already at its true sample (`rt_bids` leaves out the sync-only events). The cell counts the skipped-sample gaps (`acq_gap`) and the trials, and plots the trials. The check cell compares trial counts (and the dataset's bad channels) across all participants of the session.

**3. Breaks and cropping.** Gaps longer than 15 s between trials (surveys, results, instructions, team reassignment) are labelled `BAD_break`; anything labelled `BAD` is ignored by later steps. Then the data is cropped to the trials (motor baseline, individual, joint). The labels are added *before* cropping, while the data still start at sample 0, so their times need no correction.

**4. Bad channels.** Look at the signals and spectra, write further bad channels into `BAD_CHANNELS` in `config.py`, run the next cell. The dataset's bad channels are always included. The cell also prints the median channel SD (units check, question 18). Bad channels are **interpolated now, before ICA** (rebuilt from their neighbours), because CORRMAP needs the same channels in every participant.

**5. Save the prepared data** (500 Hz). `03_clean` applies the ICA to this file.

**6. ICA training copy.** In the reference's order: low-pass 40 Hz → downsample to 250 Hz → high-pass 1 Hz. *Why:* the high-pass removes slow drifts, which carry a lot of power and would waste ICA components; the copy is only for learning the unmixing.

**7. Dummy epochs and rejection.** The copy is cut into 1 s pieces (skipping breaks). Pieces with very unusual values (joint probability more than 3 SD from the mean) are dropped, so large one-off artifacts don't dominate the ICA.

**8. ICA.** Extended Infomax (= EEGLAB `runica` with `'extended', 1`). *Concept:* each electrode records a mixture of sources (brain, eyes, heart, muscles); ICA unmixes them into independent components, each with a scalp map and a time course. The number of components equals the data rank (27 minus interpolated channels); the reference's PCA to 50 dimensions is not possible with 27 channels.

**9. Save** the ICA, its component topographies as PNG (for choosing templates), and the numbers for the Methods (`_ica_info.csv`).

---

## `02_corrmap.ipynb`: artifact components in all participants

**1.** Loads the ICA of every participant.

**2. Templates.** You choose, by hand, one clear example per artifact type from the PNGs and write it into `CORRMAP_TEMPLATES` in `config.py`:
- **Eye blink:** strong over the forehead, both sides the same sign.
- **Lateral eye movement:** front left vs. front right, opposite signs.
- **Heartbeat:** broad diagonal gradient; regular spikes about once per second.
- **Muscle** (added to the reference): small focal spot at the edge of the head (temples, neck, jaw), often one-sided; power stays high above ~20 Hz. Add a second template if one side is missed.

MNE counts components from **0** (EEGLAB from 1).

**3. CORRMAP.** For each template, finds components with a similar scalp map (correlation ≥ 0.8) in every participant, in two rounds (template, then the average of the first-round matches).

**4. Save** per participant: `_badcomps.csv` (CORRMAP matches plus `EXTRA_EXCLUDE` from `config.py`) and `_badcomps.png`.

---

## `03_clean.ipynb`: cleaning per participant

**1.** Loads the prepared data, the ICA and the components to remove.

**2. Remove the components** from the 500 Hz data. The ICA was learned on the 250 Hz copy; its weights are a fixed spatial unmixing, so they apply to the original data (as the reference writes the weights back to the raw data). The overlay plot shows before vs. after: blinks should disappear, the rest barely change.

**3. Filter** as in `ana03`: low-pass 40 Hz (order 100), then high-pass 0.1 Hz (order 500). No re-referencing: the data stay referenced to the mastoids (the recording reference). (The 0.1 Hz high-pass barely filters with this order; see question 9.)

**4. Save** the clean data (`_clean_raw.fif`, input for the analysis), an HTML report (data and spectra before/after, ICA), and a row in `preprocessing_summary.csv` (bad channels, rejected pieces, components, duration) for the Methods.

**Next (analysis notebook):** alpha and beta envelopes on the continuous data, cut into trials relative to the same event in every participant's file, teammate/opponent labels per trial from the `trials` table, then cRQA. The reference's ERP epoching is not used.

---

## Order of work

1. Run `01_ica` for every participant (fill in bad channels as you go).
2. Look at the component PNGs, choose the templates, run `02_corrmap`. Review `_badcomps.png`; add anything missed to `EXTRA_EXCLUDE` and re-run `02_corrmap`.
3. Run `03_clean` for every participant.
