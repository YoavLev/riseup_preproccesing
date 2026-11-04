"""Python versions of the two EEGLAB functions used in the reference pipeline
(Stropahl et al., 2018, supplementary material): pop_firws (filtering) and pop_jointprob
(rejecting unusual epochs before ICA).
"""

import numpy as np
from scipy.signal import firwin, oaconvolve


def firws_filter(raw, cutoff, kind, order):
    """Filter like EEGLAB's pop_firws: windowed-sinc FIR filter with a Hann window.

    cutoff: the cut-off frequency in Hz (the -6 dB point).
    kind:   "lowpass" or "highpass".
    order:  filter order (even number); the filter has order + 1 coefficients.

    The filter is applied once and is zero-phase (its delay is compensated), so events stay
    at the same time. EEG and EOG channels are filtered in place.
    """
    if order % 2:
        raise ValueError("the filter order must be even")
    coefficients = firwin(order + 1, cutoff, window="hann", pass_zero=(kind == "lowpass"),
                          fs=raw.info["sfreq"])
    half = order // 2

    def apply(signal):
        # Repeat the first/last sample at the edges, so the edges are not pulled towards zero
        padded = np.pad(signal, half, mode="edge")
        return oaconvolve(padded, coefficients, mode="valid")

    raw.apply_function(apply, picks=["eeg", "eog"], channel_wise=True)
    return raw


def _sample_probability(x, bins=1000):
    """For every sample, the probability of its value: a histogram with `bins` bins over the
    range of the data (EEGLAB's realproba)."""
    lowest, highest = x.min(), x.max()
    index = np.floor((x - lowest) / (highest - lowest) * (bins - 1)).astype(int)
    counts = np.bincount(index, minlength=bins)
    return counts[index] / x.size


def _zscore(values):
    """Standardize across epochs (last axis)."""
    return (values - values.mean(axis=-1, keepdims=True)) / values.std(axis=-1, ddof=1, keepdims=True)


def jointprob_keep(epochs_data, threshold):
    """Which epochs to keep, like EEGLAB's pop_jointprob (local and global, threshold in SD).

    epochs_data: array (n_epochs, n_channels, n_times).

    Idea: values that rarely occur in a channel have a low probability. An epoch full of
    improbable values (a large artifact) gets a high "joint probability" score, -sum(log p).
    - local:  the score per channel, standardized across epochs;
    - global: the same score for all channels together.
    An epoch is rejected if its local score in any channel, or its global score, is more
    than `threshold` standard deviations from the mean.
    """
    n_epochs, n_channels, n_times = epochs_data.shape

    local = np.zeros((n_channels, n_epochs))
    for channel in range(n_channels):
        p = _sample_probability(epochs_data[:, channel, :].ravel())
        local[channel] = -np.log(p).reshape(n_epochs, n_times).sum(axis=1)
    local = _zscore(local)

    p_all = _sample_probability(epochs_data.ravel())
    global_ = _zscore(-np.log(p_all).reshape(n_epochs, n_channels * n_times).sum(axis=1))

    reject = (np.abs(local) > threshold).any(axis=0) | (np.abs(global_) > threshold)
    return ~reject
