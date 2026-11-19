"""Cross-recurrence quantification analysis (cRQA) of two time series, plus the tools to choose its parameters.

Used by 05_crqa.ipynb on the amplitude envelopes of two participants. Written out in numpy (no RQA package),
so every step can be explained in the Methods. References for the method: Marwan et al. (2007), Coco & Dale
(2014) (unverified, check before citing); Coco et al. (2021), the crqa R package (notes/papers/coco2021_crqa.md).

The idea in four steps:
1. **z-score** each series (mean 0, SD 1), so the radius means the same for both people.
2. **Embed** each series: each time point t becomes a vector (x[t], x[t+tau], ..., x[t+(m-1)tau]) of m values,
   tau samples apart. This turns a series into a trajectory, so "the same state" means the same short pattern,
   not just the same value. tau comes from the auto mutual information (AMI), m from false nearest neighbours (FNN).
3. **Cross-recurrence:** point (i, j) is recurrent when person A's state at time i is within the radius of
   person B's state at time j (Euclidean distance).
4. **Measures:** %REC = share of recurrent points; %DET = share of recurrent points that lie on diagonal lines of
   at least `min_line` points (both people going through the same sequence of states, possibly with a lag).
   In 05_crqa both are taken only near lag 0 (a band of diagonals), and the diagonal profile is saved.

The full N x N matrix (N = 7500 for a 300 s trial at 25 Hz) is never stored: it is computed one diagonal at a time.
"""

import numpy as np
from scipy.spatial import cKDTree


def zscore(x):
    """Mean 0, SD 1."""
    x = np.asarray(x, dtype=float)
    return (x - x.mean()) / x.std()


def embed(x, m, tau):
    """Delay embedding: rows are (x[t], x[t+tau], ..., x[t+(m-1)tau]); shape (len(x) - (m-1)*tau, m)."""
    n = len(x) - (m - 1) * tau
    if n <= 0:
        raise ValueError(f"series of {len(x)} samples is too short for m={m}, tau={tau}")
    return np.column_stack([x[i * tau:i * tau + n] for i in range(m)])


# --- Choosing tau: auto mutual information --------------------------------------------------------------------

def auto_mutual_information(x, max_lag, bins=16):
    """Mutual information (bits) between x[t] and x[t+lag], for lag = 0..max_lag (histogram estimate)."""
    x = np.asarray(x, dtype=float)
    edges = np.histogram_bin_edges(x, bins=bins)
    codes = np.clip(np.digitize(x, edges[1:-1]), 0, bins - 1)
    mi = []
    for lag in range(max_lag + 1):
        a, b = codes[:len(codes) - lag], codes[lag:]
        joint = np.zeros((bins, bins))
        np.add.at(joint, (a, b), 1)
        joint /= joint.sum()
        pa, pb = joint.sum(axis=1), joint.sum(axis=0)
        nz = joint > 0
        mi.append(np.sum(joint[nz] * np.log2(joint[nz] / np.outer(pa, pb)[nz])))
    return np.array(mi)


def first_minimum(values):
    """Index of the first local minimum (the usual choice of tau); the global minimum if there is none."""
    for i in range(1, len(values) - 1):
        if values[i] < values[i - 1] and values[i] <= values[i + 1]:
            return i
    return int(np.argmin(values[1:]) + 1)


# --- Choosing m: false nearest neighbours ---------------------------------------------------------------------

def false_nearest_neighbours(x, tau, max_dim, rtol=15.0, atol=2.0):
    """Share of false nearest neighbours for m = 1..max_dim (Kennel et al., 1992; unverified, check before citing).

    A neighbour in m dimensions is "false" when adding the (m+1)-th coordinate moves it far away: then the two
    points were only close because the trajectory was squeezed into too few dimensions.
    """
    x = zscore(x)
    shares = []
    for m in range(1, max_dim + 1):
        big = embed(x, m + 1, tau)                 # m+1 coordinates; the first m are the m-dimensional embedding
        small = big[:, :m]
        dist, idx = cKDTree(small).query(small, k=2)
        dist, idx = dist[:, 1], idx[:, 1]          # nearest neighbour other than the point itself
        extra = np.abs(big[:, m] - big[idx, m])    # how far apart the new coordinate puts them
        ok = dist > 0
        false = (extra[ok] / dist[ok] > rtol) | (np.sqrt(dist[ok] ** 2 + extra[ok] ** 2) / x.std() > atol)
        shares.append(false.mean())
    return np.array(shares)


def choose_dimension(fnn, threshold):
    """Smallest m whose share of false neighbours is at or below the threshold; else the m with the lowest share."""
    below = np.where(fnn <= threshold)[0]
    return int(below[0] + 1) if below.size else int(np.argmin(fnn) + 1)


# --- Choosing the radius --------------------------------------------------------------------------------------

def radius_for_rec(x, y, m, tau, target_rec, n_pairs=200_000, seed=0):
    """Radius at which about `target_rec` % of the cross-recurrence points are recurrent (random sample of points)."""
    a, b = embed(zscore(x), m, tau), embed(zscore(y), m, tau)
    n = min(len(a), len(b))
    rng = np.random.default_rng(seed)
    i, j = rng.integers(0, n, n_pairs), rng.integers(0, n, n_pairs)
    return float(np.percentile(np.linalg.norm(a[i] - b[j], axis=1), target_rec))


# --- Cross-recurrence measures --------------------------------------------------------------------------------

def _diagonal(a, b, k):
    """Distances on diagonal k of the cross-recurrence plot: A at time i against B at time i + k (k = lag)."""
    n = len(a)
    if k >= 0:
        return np.linalg.norm(a[:n - k] - b[k:], axis=1)
    return np.linalg.norm(a[-k:] - b[:n + k], axis=1)


def crqa(x, y, m, tau, radius, min_line=2, max_lag=None):
    """%REC and %DET (plus line statistics) of the cross-recurrence of x and y.

    x and y must have the same length (the same window of both people). Both are z-scored and embedded here.
    max_lag: only the diagonals with |lag| <= max_lag samples are used (the band around the line of
    synchronisation, lag 0); None = the whole plot. With long, smooth envelopes the whole plot cannot tell
    coupled from unrelated signals: coupling shows only near lag 0 (Coco et al., 2021: diagonal profile).
    min_line: shortest diagonal line counted for %DET, in samples (choose it in seconds: smooth envelopes make
    short lines even for unrelated signals).
    Returns a dict: rec (%), det (%), mean_line, max_line (samples), n_recurrent, n_points.
    """
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    a, b = embed(zscore(x), m, tau), embed(zscore(y), m, tau)
    n = len(a)
    lags = range(-(n - 1), n) if max_lag is None else range(-max_lag, max_lag + 1)
    n_recurrent = n_points = 0
    lines = []                                     # lengths of all diagonal lines >= min_line
    for k in lags:
        rec = _diagonal(a, b, k) <= radius
        n_points += rec.size
        count = int(rec.sum())
        if count == 0:
            continue
        n_recurrent += count
        # lengths of the runs of consecutive recurrent points on this diagonal
        edges = np.diff(np.concatenate([[0], rec.astype(np.int8), [0]]))
        lengths = np.where(edges == -1)[0] - np.where(edges == 1)[0]
        lines.extend(lengths[lengths >= min_line])
    lines = np.array(lines)
    on_lines = int(lines.sum()) if lines.size else 0
    return {
        "rec": 100 * n_recurrent / n_points,
        "det": 100 * on_lines / n_recurrent if n_recurrent else np.nan,
        "mean_line": float(lines.mean()) if lines.size else np.nan,
        "max_line": int(lines.max()) if lines.size else 0,
        "n_recurrent": n_recurrent,
        "n_points": n_points,
    }


def recurrence_profile(x, y, m, tau, radius, max_lag):
    """Diagonal recurrence profile: %REC on each diagonal, for lags -max_lag..+max_lag samples.

    Positive lag k: A at time i matches B at time i + k, i.e. B follows A by k samples.
    Returns (lags, rec_percent).
    """
    a, b = embed(zscore(x), m, tau), embed(zscore(y), m, tau)
    lags = np.arange(-max_lag, max_lag + 1)
    return lags, np.array([100 * np.mean(_diagonal(a, b, k) <= radius) for k in lags])
