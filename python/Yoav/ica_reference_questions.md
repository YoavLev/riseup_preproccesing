# Open questions: adapting the reference ICA pipeline

**Reference:** Stropahl, Bauer, Debener & Bleichner (2018), supplementary material (`notes/papers/eeg_data_sheet.PDF`), scripts `ana01_ICA.m`, `ana02_corrmap.m`, `ana03_preprocessing.m`. The notebooks `01_ica`, `02_corrmap`, `03_clean` follow it: same order, same constants (in `config.py`).

Below: where our data don't fit the reference. The numbers are kept from earlier versions because `notes/decisions.md` refers to them.

---

## Open: needs a decision from the supervisor

### 1. Filters of the ICA copy: 1–40 Hz
I was asked to "change the high-pass to 40 Hz". In the reference, 40 Hz is the **low-pass**; a 40 Hz high-pass would remove alpha and beta. I used 1–40 Hz, as in the reference. **Confirm.**

### 5. Bad channels interpolated before ICA
CORRMAP needs the same channels in every participant, so bad channels are rebuilt from their neighbours **before** ICA. The reference has no bad-channel step. (Session01: no bad channels.) **OK?**

### 6. No PCA reduction, no average reference
The reference reduces the data to 50 dimensions before ICA. We have only 27 EEG channels, so the code caps it at the data rank: 27 components, i.e. **no reduction**. Each interpolated channel lowers the rank by 1. We don't average-reference (item 8). If we did, the rank would drop by 1 more and the ICA would need one component less. **OK, or reduce further?**

### 7. EOG channels left out of the ICA
The reference decomposes all channels. Our EOG2–4 have no scalp positions, so ICA uses the 27 EEG channels only. Eye components are still found by their maps. **Include the EOG channels?**

### 9. Filter lengths were set for 1000 Hz data; ours are 500 Hz
A filter order is counted in samples, so the same number gives a different filter here. The constants are kept as written. Effects:
- 40 Hz low-pass: starts cutting a bit later (~32 Hz instead of ~25 Hz). Beta (≤ 30 Hz) is untouched.
- **0.1 Hz high-pass: barely filters** (it keeps 90% of the offset at 0 Hz, also in the reference). The filter is too short for such a low frequency.
- **Does it matter?** Not for the envelopes: they are band-passed at 8–12 / 13–30 Hz, which removes the offset. It only makes the "segment image" panel in the component PNGs useless.

**Keep as is, or make the high-pass work (longer filter, or remove the mean first)?**

### 10. Rejecting unusual 1 s pieces before ICA (Python version)
MNE has no `pop_jointprob`, so `utils.jointprob_keep` rebuilds EEGLAB's method (3 SD). Numbers can differ slightly from EEGLAB. Session01 rejected 3–16% of the pieces (172–778 of 4964). **Acceptable?**

### 14. No extra artifact rejection after cleaning
The reference cuts short ERP epochs and rejects them at 4 SD. We keep whole trials of continuous data (cRQA needs long segments), so this step is not used. **Should leftover artifacts be marked as bad segments?**

### 15. CORRMAP templates for new sessions
The session01 templates are **provisional** (suggested by ICLabel): blink sub-005 #4, lateral eye sub-003 #10, muscle sub-004 #21. **No heartbeat template**: no component looked clearly like heartbeat (best ICLabel probability 0.49). **OK to use ICLabel for suggestions and to leave out heartbeat?** For sessions 02–03, rerun CORRMAP over all 18 participants **with the same templates**, so earlier choices stay (mostly) the same.

### 17. Alpha peak
In the pilot, parietal/occipital channels showed no clear alpha peak. That can be normal in an active visual game, but alpha is one of our two bands. **Check in session01** (component spectra in `04`, motor baseline and breaks).

### 20. Automatic bad-channel criteria (new)
`utils.find_bad_channels` uses PREP's criteria, adapted to our 27-channel cap: noisy = robust amplitude z > 5; flat = < 1 µV; uncorrelated = highest correlation with any other channel < 0.4 in **> 50%** of 1 s windows (PREP: 1%, which flagged 2–11 channels per person here). Fp1/Fp2 are not checked for noise (their large signal is eye activity, which ICA removes). **OK?**

---

## Settled (no action needed)

- **2.** Mastoids M1/M2 removed (supervisor).
- **3.** Motor baseline kept in the data.
- **4.** Muscle template added to CORRMAP (decided 2026-10-29).
- **8.** Reference: the data stay referenced to the mastoids (the recording reference); no re-referencing (decided 2026-10-29).
- **11.** Breaks between trials (> 15 s) are skipped when cutting the 1 s pieces for ICA.
- **12.** MNE's CORRMAP has no "at most 3 components per participant" limit. Removed per participant in session01: 2–4.
- **13.** Component numbers start at 0 in MNE (1 in EEGLAB).
- **18.** Units are correct in the BIDS data: median channel SD 7–12 µV (closed 2026-11-05).
- **19.** Channels marked bad in the dataset's `channels.tsv` are always treated as bad.
