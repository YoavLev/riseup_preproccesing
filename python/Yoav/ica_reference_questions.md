# Open questions: adapting the reference ICA pipeline

**Reference:** Stropahl, Bauer, Debener & Bleichner (2018), *Source-modelling auditory processes of EEG data using EEGLAB and Brainstorm*, supplementary material (`notes/papers/eeg_data_sheet.PDF`), scripts `ana01_ICA.m`, `ana02_corrmap.m`, `ana03_preprocessing.m`.

**What was done:** the preprocessing was rebuilt to follow that pipeline: same order, same constants, same structure (three steps → three notebooks: `01_ica`, `02_corrmap`, `03_clean`; all constants in `config.py`).

Some things in the reference don't transfer directly to the RiseTogether data. Each one is listed below: **the issue**, **what I did for now**, and **the question to decide**. Items 1–3 were requested; items 4 onward need a decision or a confirmation.

---

## Requested changes

### 1. "Change the high-pass to 40 Hz"
- **Issue:** in the reference, 40 Hz is the **low-pass** for the ICA copy (`LP = 40`); the high-pass is 1 Hz (`HP = 1`). A *high*-pass at 40 Hz would remove everything below 40 Hz, including alpha and beta, so it can't be what was meant.
- **What I did:** the ICA copy is filtered **1–40 Hz**, as in the reference (low-pass 40 Hz, high-pass 1 Hz).
- **Question:** confirm this is what was meant.

### 2. Mastoids removed
- **What I did:** M1 and M2 are removed right after loading (`DROP_CHANNELS` in `config.py`). The data has 27 scalp EEG channels + 3 EOG channels.
- **Consequence:** the data can no longer be re-referenced to the mastoids. See item 8 (reference).

### 3. Motor baseline kept
- **What I did:** the motor-baseline trial is kept in the data (`CONDITIONS` in `config.py`), with the individual and joint trials.

---

## Needs a decision or confirmation

### 4. Muscle artifacts: no template, and the 40 Hz low-pass hides them — **decided 2026-10-29: muscle template added**
- **Issue:** muscle activity (EMG) is mostly above 20 Hz and extends well above 40 Hz. With the ICA copy low-passed at 40 Hz, ICA sees less of it, and the reference has CORRMAP templates only for **blinks, heartbeat and lateral eye movements**, not muscle. In this task, participants tap buttons the whole time, and **beta (13–30 Hz) is one of our two bands**, so leftover muscle activity matters more here than in the reference (an auditory listening task).
- **What I did:** followed the reference (3 templates, 40 Hz). `EXTRA_EXCLUDE` in `config.py` lets you remove extra components after looking at them.
- **Decision:** a 4th CORRMAP template, `muscle`, was added (`CORRMAP_TEMPLATES` in `config.py`). Muscle components are focal and often one-sided, so a second template (left/right) may be needed; check `_badcomps.png` after running CORRMAP.

### 5. CORRMAP needs the same channels in every participant → bad channels interpolated *before* ICA
- **Issue:** CORRMAP compares component maps *across* participants, so every participant must have exactly the same channels (MNE refuses otherwise; this was tested). The reference has no bad-channel step at all.
- **What I did:** bad channels are interpolated (rebuilt from their neighbours) **before** ICA, and the number of ICA components is reduced by the number of interpolated channels (they add no information).
- **Question:** OK? (The alternative, interpolating after ICA, is incompatible with CORRMAP.)

### 6. PCA with 50 dimensions is impossible with 27 channels
- **Issue:** the reference reduces the data to 50 dimensions with PCA before ICA (`PCA = 1`, `PCADIMS = 50`); their cap had more channels. We have 27 EEG channels, so 50 dimensions don't exist.
- **What I did:** kept the setting, capped at the data rank: 27 components (fewer if channels were interpolated). In practice: **no PCA reduction**.
- **Question:** OK, or reduce further (e.g., 20 components)?

### 7. EOG channels are not part of the ICA
- **Issue:** the reference runs ICA on all recorded channels. We have 3 EOG channels (EOG2–4) without standard scalp positions.
- **What I did:** ICA uses the 27 EEG channels only (EOG channels are kept in the data but not decomposed). Eye artifacts are still found through their scalp maps (CORRMAP).
- **Question:** include the EOG channels in the ICA?

### 8. Reference: the data stay in the recording reference — **decided 2026-10-29: recording reference = mastoids, kept**
- **Issue:** the reference pipeline never re-references. With the mastoids removed, the options are the **recording reference** (whatever the hyperscanners used during recording) or the **average reference**. We don't yet know what the recording reference of the hyperscanners is.
- **What I did:** no re-referencing, as in the reference pipeline.
- **Decision:** the recording reference is the **mastoids** (confirmed). The data stay in that reference (no re-referencing), as in the reference pipeline; the M1/M2 channels themselves are removed from the data.

### 9. The filter constants were chosen for 1000 Hz data; ours is 500 Hz
- **Issue:** a filter "order" is counted in samples, so the same order gives a different filter at a different sampling rate (the reference itself notes "filter order depends on sampling rate"). Their data: 1000 Hz. Ours: 500 Hz. I kept the constants exactly as written. Measured effect:

| Filter | In the reference (1000 Hz) | Here (500 Hz) |
|---|---|---|
| ICA low-pass 40 Hz, order 100 | starts cutting at ~25 Hz (−0.9 dB at 30 Hz) | starts cutting at ~32 Hz (0 dB at 30 Hz): steeper |
| ICA high-pass 1 Hz, order 500 (applied at 250 Hz in both) | identical | identical |
| Analysis low-pass 40 Hz, order 100 | as first row | as first row; beta (up to 30 Hz) untouched |
| Analysis high-pass 0.1 Hz, order 500 | **barely filters**: −0.4 dB even at 0 Hz | **barely filters**: −0.9 dB at 0 Hz |

- **The 0.1 Hz high-pass is effectively not a filter**, in the reference and here: the filter is far too short for such a low cut-off, so slow drifts pass almost unchanged. For our alpha/beta envelopes this doesn't matter (the analysis band-passes 8–12 and 13–30 Hz anyway).
- **Question:** keep the constants literally (current), or adjust the orders to 500 Hz / make the 0.1 Hz high-pass real (much longer filter)?

### 10. Joint-probability rejection rewritten in Python
- **Issue:** the reference uses EEGLAB's `pop_jointprob` to throw away unusual 1 s pieces before ICA (3 SD). MNE has no equivalent.
- **What I did:** a Python version (`jointprob_keep` in `utils.py`) following EEGLAB's method: histogram-based probability of each value, per channel ("local") and over all channels ("global"), standardized across epochs, rejected above 3 SD. It is a re-implementation, not EEGLAB's code: numbers can differ slightly.
- **Question:** acceptable? (The share of rejected pieces is reported per participant.)

### 11. Dummy epochs skip the breaks
- **Issue:** the reference cuts the *whole* recording into 1 s pieces. Our recordings include surveys, instructions and team reassignment between trials.
- **What I did:** the breaks between trials (> 15 s) are skipped, so ICA learns from task data only (the joint-probability step would catch some of it anyway).

### 12. CORRMAP's "up to 3 components per participant" has no MNE equivalent
- **Issue:** the reference calls CORRMAP with `'ics', 3`: at most 3 components per participant can match a template. MNE's CORRMAP takes all components above the threshold (0.8).
- **What I did:** MNE's behaviour; the number of removed components per participant is reported, so an unusual count is visible.

### 13. Component numbers start at 0, not 1
- **Not a question, a warning:** EEGLAB counts components from 1, MNE from 0. The reference's example "component 1 = blink" would be component 0 here. Templates in `config.py` use MNE numbering.

### 14. The reference's epoching (ana03) is not used
- **Issue:** `ana03` cuts −0.2 to 0.8 s epochs around a stimulus and rejects them at 4 SD. This thesis analyzes **whole trials** of continuous data with cRQA, which needs long segments.
- **What I did:** filtering of `ana03` kept; epoching and epoch rejection not used. Trials are cut in the analysis notebook.
- **Question:** is another artifact-rejection step wanted for the continuous data (e.g., marking residual artifacts as bad segments)?

### 15. Templates and new sessions
- **Issue:** CORRMAP picks templates from one participant and searches all others. The current pilot has 5 participants with EEG (Padme's headset could not be aligned).
- **Question:** run CORRMAP once over all sessions at the end of data collection, or per session? If it is re-run when sessions are added, keep the same templates, so the components already chosen for earlier participants stay the same. (CORRMAP's second round uses the average of all matches, so adding participants can still shift borderline cases slightly.)

---

## Observations from the first test run (pilot session `SCOOL`, 5 participants)

The three notebooks were run end to end on the pilot data, with **placeholder** CORRMAP templates (so which components were removed means nothing yet).

### 16. A large share of 1 s pieces is rejected before ICA
- Rejected by the joint-probability step (3 SD): Parsnip 28%, Patrick 9%, Peter 16%, Pia 31%, Polly 16%.
- Not a problem for ICA (700–1000 s of data remain per participant), but higher than one might expect. Possibly a lot of movement during play, or a sign that the Python version of `pop_jointprob` (item 10) is stricter than EEGLAB's.

### 17. No clear alpha peak in parietal/occipital channels, even before ICA
- In all 5 participants, power in 8–12 Hz is *lower* than in the neighbouring 5–7 and 14–17 Hz, before and after cleaning (so ICA did not remove alpha).
- This can be normal during an active, eyes-open visual game (alpha is suppressed while attending), but it matters for RQ1, which relies on alpha. Worth checking the motor baseline and the breaks (where alpha may be stronger), and the recording reference (item 8).

### 18. The stored EEG values are about 1000× too large (units)
- **Issue:** in all 5 files exported by `build.py`, a typical channel varies by about 7,000–13,000 "µV" (robust standard deviation, 1–40 Hz). Real scalp EEG varies by about 5–50 µV. The values look like **nanovolts stored as if they were microvolts** (`eeg_export.py` multiplies by `1e-6`, i.e. assumes the hyperscanner writes µV).
- **Effect on this pipeline:** none on the cleaning itself (ICA, CORRMAP and the joint-probability step don't depend on the absolute scale), and none on cRQA if envelopes are standardized. But plots, the MNE report and any amplitude reported in µV are 1000× off.
- **Question for Luke:** in which unit does the hyperscanner write EEG? If nV, the export scale should be `1e-9`. (Not changed here: it's Luke's code, and it affects every student using the files.)
- **Update 2026-11-04 (BIDS dataset):** the BDF files store the ADC counts × the ADC step in µV, and Luke's example alpha PSDs (~1e-12 V²/Hz) are a plausible size, so this is probably fixed in the BIDS export. **To verify once** in `01_ica` step 4, which prints the median channel SD (1–40 Hz): it should be tens of µV. Then close this item.
- Also seen: one channel in Parsnip's file is flat (SD 0) and one in Pia's is ~6× noisier than the rest: likely bad channels for step 4 of `01_ica`.

---

## New with the BIDS dataset (2026-11-04)

### 19. Bad channels from the dataset
- Channels at the ADC limit more than 50% of the time are `status = bad` in `channels.tsv`. `01_ica` always treats them as bad (and interpolates them), plus the ones in `BAD_CHANNELS`. Both lists are in `_ica_info.csv` (`bad_channels_dataset`, `bad_channels`).
