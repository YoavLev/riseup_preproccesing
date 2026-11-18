# Open questions: adapting the reference ICA pipeline

**Reference:** Stropahl, Bauer, Debener & Bleichner (2018), supplementary material (`notes/papers/eeg_data_sheet.PDF`), scripts `ana01_ICA.m`, `ana02_corrmap.m`, `ana03_preprocessing.m`. The notebooks `01_ica`, `02_corrmap`, `03_clean` follow it: same order, same constants (in `config.py`).

Below: where our data don't fit the reference. The numbers stay fixed because `notes/decisions.md` refers to them.

---

## Open: needs a decision from the supervisor

### 1. Filters of the ICA copy: 1–40 Hz
I was asked to "change the high-pass to 40 Hz". In the reference, 40 Hz is the **low-pass**; a 40 Hz high-pass would remove alpha and beta. I used 1–40 Hz, as in the reference. **Confirm.**

### 5. Bad channels interpolated before ICA
CORRMAP needs the same channels in every participant, so bad channels are rebuilt from their neighbours **before** ICA. The reference has no bad-channel step. **OK?**

### 6. No PCA reduction
The reference reduces the data to 50 dimensions before ICA. With 27 EEG channels the code caps it at the data rank: 27 − interpolated channels − 1 (common average reference), i.e. **no real reduction** (26 components without bad channels). **OK, or reduce further?**

### 7. EOG channels left out of the ICA
The reference decomposes all channels. Our EOG2–4 have no scalp positions, so ICA uses the 27 EEG channels only. Eye components are still found by their maps. **Include the EOG channels?**

### 9. MNE's `filter()` instead of the reference's fixed-order filters (changed 2026-11-12, confirm)
Same cut-offs, but MNE chooses the filter length. The recordings drift slowly by several mV (sub-016: C4 6.8 mV SD before filtering):
- The reference's 1 Hz high-pass (order 500 at 250 Hz, 2 s long) lets ~1/80 of the drift through (−38 dB): ~80 µV on C4, more than the EEG; about ten of 26 ICA components were copies of it.
- The 0.1 Hz high-pass (order 500 at 500 Hz) attenuates only 0.9 dB below 0.5 Hz: it did not filter.
- With MNE's filter (75–300 dB near 0 Hz), sub-016's ICA had 15 ICLabel brain components instead of 5 (`scripts/ica_test.py sub-016 mne-filter`).

### 10. Rejecting unusual 1 s pieces before ICA (Python version)
MNE has no `pop_jointprob`, so `utils.jointprob_keep` rebuilds EEGLAB's method (3 SD). Numbers can differ slightly from EEGLAB. Rejected: 3–16% of the pieces (session01). **Acceptable?**

### 14. Artifact rejection after cleaning
The reference cuts short ERP epochs and rejects them at 4 SD. We keep whole trials of continuous data (cRQA needs long segments). Stretches where all channels jump are marked `BAD_amplitude` in `01` (step 3b, 1 s windows, median channel SD > 50 µV). **Open: when is a whole trial dropped (share of bad time)?**

### 15. CORRMAP templates per session
Templates are chosen per session from that session's participants, suggested by ICLabel (`scripts/suggest_corrmap_templates.py`) and checked by eye. Session03: blink sub-013 #1, lateral eye sub-014 #1. **No heartbeat template** (no clear heart component). **OK to use ICLabel for suggestions, to leave out heartbeat, and to choose templates per session** (rather than one set for all sessions)?

### 17. Alpha peak
In the pilot, parietal/occipital channels showed no clear alpha peak. That can be normal in an active visual game, but alpha is one of our two bands. **Check** in the spectra of the session and localizer components (`04`, `04a`).

### 20. Automatic bad-channel criteria
`utils.find_bad_channels` uses PREP's criteria, adapted to our 27-channel cap: noisy = robust amplitude z > 5; flat = < 1 µV; uncorrelated = highest correlation with any other channel < 0.4 in **> 50%** of 1 s windows (PREP: 1%, which flagged 2–11 channels per person here). Fp1/Fp2 are not checked for noise (their large signal is eye activity, which ICA removes). **OK?**

### 21. ICA does not converge within 500 iterations
MNE's extended Infomax stops at 500 iterations (EEGLAB: 512). Nearly all session and localizer ICAs reach the cap. **Raise the cap (e.g. 2000) and check whether the maps change?**

---

## Settled (no action needed)

- **2.** Mastoids M1/M2 removed (supervisor).
- **3.** Motor baseline kept in the data.
- **4.** No muscle template (supervisor, 2026-11-11: removing muscle components can also remove brain signal). Muscle stays in the data: watch the beta band.
- **8.** Common average reference, applied in `01` after bad-channel interpolation (supervisor, decided 2026-11-11). The recording reference was the mastoids.
- **11.** Breaks between trials (> 15 s) are skipped when cutting the 1 s pieces for ICA.
- **12.** MNE's CORRMAP has no "at most 3 components per participant" limit.
- **13.** Component numbers start at 0 in MNE (1 in EEGLAB).
- **18.** Units are correct in the BIDS data: median channel SD 7–12 µV (closed 2026-11-05).
- **19.** Channels marked bad in the dataset's `channels.tsv` are always treated as bad.
