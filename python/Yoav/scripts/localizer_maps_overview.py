"""One overview PNG per localizer: the scalp maps of every component of every participant of a session.

Rows = participants, columns = components (ICA000, ICA001, ...) of that participant's localizer ICA
(04a_localizer_ica). Each map has its own colour scale (a component's sign and size are arbitrary).
Makes it easy to compare the participants and to find a similar visual or motor component in each.

Output: <OUT_ROOT>/spatial_filter_review/<session>/localizer/<session>_<localizer>_all_maps.png

Usage (on UCloud, in the riseup_preproccesing folder):
    uv run python python/Yoav/scripts/localizer_maps_overview.py session03
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # no screen needed
import matplotlib.pyplot as plt
import mne

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # python/Yoav
import config as cfg  # noqa: E402


def overview(session, localizer):
    """Draw and save the overview of one localizer; returns the path of the PNG."""
    subjects = cfg.session_subjects(session)
    icas = {}
    for sub in subjects:
        path = cfg.OUT_ROOT / sub / f"{sub}_{localizer}_localizer_ica.fif"
        if path.exists():
            icas[sub] = mne.preprocessing.read_ica(path, verbose="error")
        else:
            print(f"{sub}: no {localizer} localizer ICA (run 04a_localizer_ica first)")
    if not icas:
        return None

    n_cols = max(ica.n_components_ for ica in icas.values())
    fig, axes = plt.subplots(len(icas), n_cols, figsize=(1.1 * n_cols, 1.3 * len(icas)), squeeze=False)
    for row, (sub, ica) in enumerate(icas.items()):
        maps = ica.get_components()  # channels x components (mixing matrix)
        for comp in range(n_cols):
            ax = axes[row, comp]
            if comp >= ica.n_components_:  # fewer components (interpolated channels): empty cell
                ax.axis("off")
                continue
            mne.viz.plot_topomap(maps[:, comp], ica.info, axes=ax, show=False, contours=0, sensors=False)
            ax.set_title(f"{comp:03d}", fontsize=7, pad=1)
        axes[row, 0].text(-0.35, 0.5, sub, transform=axes[row, 0].transAxes, rotation=90,
                          ha="center", va="center", fontsize=9)

    fig.suptitle(f"{session}: ICA on the {localizer} localizer, all components (rows: participants)", fontsize=11)
    fig.tight_layout(rect=(0.01, 0, 1, 0.97))
    out = cfg.OUT_ROOT / "spatial_filter_review" / session / "localizer" / f"{session}_{localizer}_all_maps.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=110)
    plt.close(fig)
    return out


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    for name in cfg.LOCALIZERS:
        print(name, "->", overview(sys.argv[1], name))
