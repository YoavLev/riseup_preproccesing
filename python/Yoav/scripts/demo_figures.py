"""Figures for the cRQA demo (2026-11-24): the steps from EEG to envelope, for the participants in SPATIAL_FILTER_CHOICE.

Saves to <derivatives>/demo/<session>/:
  eeg_sample_<sub>.png        10 s of the clean EEG (all scalp channels) during a Team trial
  components_<session>.png    scalp map + spectrum of each chosen component (rows = participants)
  envelope_steps_<sub>_<area>_<seconds>s.png   as 04 section 3b: component -> band-pass -> Hilbert envelope ->
                              saved 25 Hz envelope -> smoothed at 1 Hz (zero-phase FIR, Hamming, as in 05)

Usage (on UCloud, in riseup_preproccesing):
    uv run python python/Yoav/scripts/demo_figures.py session03
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mne
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as cfg  # noqa: E402

mne.set_log_level("error")
SESSION = sys.argv[1] if len(sys.argv) > 1 else "session03"
OUT = cfg.OUT_ROOT.parent / "demo" / SESSION
OUT.mkdir(parents=True, exist_ok=True)
SMOOTH_HZ = 1.0
EEG_SECONDS = 10        # length of the EEG sample
EEG_OFFSET = 60         # s after the start of the first Team trial
STEP_T0 = 450           # start of the envelope-step plots [s from the start of the clean data], as in 04 3b
subjects = [s for s in cfg.session_subjects(SESSION) if s in cfg.SPATIAL_FILTER_CHOICE]


def load(sub):
    raw = mne.io.read_raw(cfg.out_dir(sub) / f"{sub}_clean_raw.fif", preload=True)
    return raw


def ica_of(sub, name):
    path = cfg.out_dir(sub) / (f"{sub}_ica.fif" if name == "session" else f"{sub}_{name}_localizer_ica.fif")
    return mne.preprocessing.read_ica(path)


def first_joint_onset(raw):
    """Time (s from the start of the clean data) of the first Team trial start."""
    for onset, desc in zip(raw.annotations.onset, raw.annotations.description):
        if desc.startswith("TRIAL_START/joint/"):
            return onset - raw.first_time if raw.annotations.orig_time is not None else onset
    return 0.0


# --- 1. EEG sample -------------------------------------------------------------------------------------------
for sub in subjects:
    raw = load(sub)
    eeg = raw.copy().pick("eeg")
    t0 = first_joint_onset(raw) + EEG_OFFSET
    sf = eeg.info["sfreq"]
    data = eeg.get_data(start=int(t0 * sf), stop=int((t0 + EEG_SECONDS) * sf)) * 1e6   # µV
    t = np.arange(data.shape[1]) / sf
    spacing = 60  # µV between channels
    fig, ax = plt.subplots(figsize=(11, 10))
    for i, name in enumerate(eeg.ch_names):
        ax.plot(t, data[i] - i * spacing, color="k", lw=0.5)
    ax.set_yticks(-np.arange(len(eeg.ch_names)) * spacing, eeg.ch_names, fontsize=7)
    ax.set(xlabel="time (s)", xlim=(0, EEG_SECONDS + 0.6),
           title=f"{sub}: clean EEG (0.1-40 Hz, average reference, eye components removed), "
                 f"{EEG_SECONDS} s of the first Team trial")
    # scale bar to the right of the traces
    ax.plot([EEG_SECONDS + 0.3] * 2, [0, -spacing], color="C3", lw=2)
    ax.text(EEG_SECONDS + 0.38, -spacing / 2, f"{spacing} µV", color="C3", ha="left", va="center", fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / f"eeg_sample_{sub}.png", dpi=130)
    plt.close(fig)

# --- 2. The chosen components: scalp map + spectrum ---------------------------------------------------------
areas = sorted({a for s in subjects for a in cfg.SPATIAL_FILTER_CHOICE[s]}, reverse=True)   # visual, motor
fig, axes = plt.subplots(len(subjects), 2 * len(areas), figsize=(4 * len(areas) * 1.6, 3.2 * len(subjects)),
                         gridspec_kw={"width_ratios": [1, 1.6] * len(areas)}, squeeze=False)
for row, sub in enumerate(subjects):
    raw = load(sub)
    for k, area in enumerate(areas):
        ica_name, comp, band = cfg.SPATIAL_FILTER_CHOICE[sub][area]
        ica = ica_of(sub, ica_name)
        ax_map, ax_psd = axes[row, 2 * k], axes[row, 2 * k + 1]
        mne.viz.plot_topomap(ica.get_components()[:, comp], ica.info, axes=ax_map, show=False, contours=0)
        ax_map.set_title(f"{sub}  {area}\n{ica_name} ICA{comp:03d}", fontsize=9)
        src = ica.get_sources(raw).pick([f"ICA{comp:03d}"])
        psd = src.compute_psd(fmin=1, fmax=40, picks="all", n_fft=int(4 * raw.info["sfreq"]))
        ax_psd.semilogy(psd.freqs, psd.get_data()[0], color="k", lw=1)
        lo, hi = cfg.SPATIAL_FILTER_BANDS[area][band]
        ax_psd.axvspan(lo, hi, color="C0", alpha=0.15)
        ax_psd.set(xlabel="frequency (Hz)", ylabel="power (a.u.)", title=f"spectrum; shaded: {band} {lo}-{hi} Hz")
        ax_psd.title.set_fontsize(9)
        ax_psd.set_yticks([])
fig.tight_layout()
fig.savefig(OUT / f"components_{SESSION}.png", dpi=130)
plt.close(fig)

# --- 3. From component to envelope (as 04 section 3b) -------------------------------------------------------
for sub in subjects:
    raw = load(sub)
    for area, (ica_name, comp, band) in cfg.SPATIAL_FILTER_CHOICE[sub].items():
        ica = ica_of(sub, ica_name)
        lo, hi = cfg.SPATIAL_FILTER_BANDS[area][band]
        step1 = ica.get_sources(raw).pick([f"ICA{comp:03d}"])
        step2 = step1.copy().filter(lo, hi, picks="all", l_trans_bandwidth=cfg.SPATIAL_FILTER_TRANSITION,
                                    h_trans_bandwidth=cfg.SPATIAL_FILTER_TRANSITION)
        step3 = step2.copy().apply_hilbert(envelope=True, picks="all")
        saved = mne.io.read_raw(cfg.out_dir(sub) / f"{sub}_{area}_{band}_{ica_name}_envelope_raw.fif", preload=True)
        smooth = saved.copy().filter(None, SMOOTH_HZ, picks="all", h_trans_bandwidth=SMOOTH_HZ / 2, method="fir",
                                     phase="zero", fir_window="hamming", fir_design="firwin")
        for seconds in (5, 30):
            def win(x):
                s = x.info["sfreq"]
                a, b = int(STEP_T0 * s), int((STEP_T0 + seconds) * s)
                return x.times[a:b], x.get_data(start=a, stop=b)[0] * 1e6
            fig, axes = plt.subplots(4, 1, figsize=(11, 8.5), sharex=True)
            t, y = win(step1)
            axes[0].plot(t, y, color="0.3", lw=0.6)
            axes[0].set_title(f"1. {ica_name} ICA{comp:03d}: component time course (all frequencies)", fontsize=10)
            t, y = win(step2)
            axes[1].plot(t, y, color="C0", lw=0.6)
            axes[1].set_title(f"2. band-passed {lo}-{hi} Hz ({band})", fontsize=10)
            axes[2].plot(t, y, color="C0", lw=0.5, alpha=0.5, label="band-passed")
            t, y = win(step3)
            axes[2].plot(t, y, color="C3", lw=1.3, label="Hilbert envelope (500 Hz)")
            t, y = win(saved)
            axes[2].plot(t, y, "o", color="k", ms=2, label="saved envelope (25 Hz)")
            axes[2].legend(loc="upper right", fontsize=8)
            axes[2].set_title("3. amplitude envelope (Hilbert), downsampled to 25 Hz", fontsize=10)
            axes[3].plot(t, y, color="k", lw=0.6, alpha=0.4, label="saved envelope (25 Hz)")
            t, y = win(smooth)
            axes[3].plot(t, y, color="C2", lw=2, label=f"smoothed, low-pass {SMOOTH_HZ:g} Hz")
            axes[3].legend(loc="upper right", fontsize=8)
            axes[3].set_title(f"4. smoothed envelope (zero-phase FIR, Hamming window, {SMOOTH_HZ:g} Hz)", fontsize=10)
            axes[3].set_xlabel("time from the start of the clean data (s)")
            for ax in axes:
                ax.set_ylabel("a.u.")
                ax.set_yticks([])
                for s in ("top", "right"):
                    ax.spines[s].set_visible(False)
            fig.suptitle(f"{sub}, {area}: from component to envelope ({seconds} s)", fontsize=11)
            fig.tight_layout()
            fig.savefig(OUT / f"envelope_steps_{sub}_{area}_{seconds}s.png", dpi=130)
            plt.close(fig)

print("saved to", OUT)
for p in sorted(OUT.iterdir()):
    print(" ", p.name)
