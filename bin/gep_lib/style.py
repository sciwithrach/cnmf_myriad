"""Plot style: UCL colours, palettes with no repeated colours, and publication-ready matplotlib settings."""
from __future__ import annotations

import itertools
from typing import Dict, Iterable, List

import matplotlib as mpl
import numpy as np
import scanpy as sc
import seaborn as sb

# --- UCL colour scheme ------------------------------------------------------------------------------------------
# Core
DARKPUR = '#361a54'
BRIGHTPUR = '#993bff'
MIDPUR = '#ba82ff'
# Supporting
LIGHTPUR = '#ddbdff'
PALEPUR = '#eedeff'
HERITAGEBLU = '#30d6ff'
# Graph colours, dark
DARKBLU = '#002ea6'
DARKORA = '#781c1c'
DARKFUS = '#9e1a54'
DARKTEA = '#005e5c'
# Graph colours, mid
BLU = '#5487ff'
ORA = '#E36C2A'
FUS = '#ED367D'
TEA = '#57b444'

# Ordered so that neighbours in the list are easy to tell apart
UCL_BASE = [BRIGHTPUR, ORA, TEA, BLU, FUS, HERITAGEBLU, DARKPUR, DARKORA, DARKTEA, DARKBLU, DARKFUS, MIDPUR, LIGHTPUR]

MM = 1 / 25.4                 # inches per millimetre
SINGLE_COLUMN = 85 * MM       # figure widths used by most journals
DOUBLE_COLUMN = 180 * MM


# --- Colour maths -----------------------------------------------------------------------------------------------
def _srgb_to_lab(rgb: np.ndarray) -> np.ndarray:
    """Convert an (n, 3) array of sRGB values in 0-1 to CIELAB (D65), for measuring colour differences."""
    rgb = np.asarray(rgb, dtype=float)
    lin = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    m = np.array([[0.4124564, 0.3575761, 0.1804375],
                  [0.2126729, 0.7151522, 0.0721750],
                  [0.0193339, 0.1191920, 0.9503041]])
    xyz = lin @ m.T / np.array([0.95047, 1.0, 1.08883])
    f = np.where(xyz > 216 / 24389, np.cbrt(xyz), (24389 / 27 * xyz + 16) / 116)
    return np.stack([116 * f[:, 1] - 16, 500 * (f[:, 0] - f[:, 1]), 200 * (f[:, 1] - f[:, 2])], axis=1)


def _mix(colour: str, other: str, amount: float) -> str:
    """Blend `colour` towards `other` by `amount` (0 = colour, 1 = other)."""
    a = np.array(mpl.colors.to_rgb(colour))
    b = np.array(mpl.colors.to_rgb(other))
    return mpl.colors.to_hex(a * (1 - amount) + b * amount)


def ucl_palette(n: int, min_distance: float = 12.0) -> List[str]:
    """Return `n` colours built from the UCL scheme, all different from each other.

    The pool is the UCL colours plus lighter and darker versions of each. Colours are picked one at a time,
    always taking the one furthest (in CIELAB) from those already chosen, so no two clusters share a colour.
    A warning is printed if the closest pair is still closer than `min_distance` (a just-noticeable difference is ~2).
    """
    pool: List[str] = list(UCL_BASE)
    for amount in (0.35, 0.6):
        pool += [_mix(c, '#ffffff', amount) for c in UCL_BASE] + [_mix(c, '#000000', amount) for c in UCL_BASE]
    pool = list(dict.fromkeys(pool))
    lab = _srgb_to_lab(np.array([mpl.colors.to_rgb(c) for c in pool]))

    chosen = [0]
    while len(chosen) < min(n, len(pool)):
        dist = np.min(np.linalg.norm(lab[:, None, :] - lab[chosen][None, :, :], axis=2), axis=1)
        dist[chosen] = -1
        chosen.append(int(np.argmax(dist)))

    colours = [pool[i] for i in chosen]
    if n > len(pool):
        raise ValueError(f'Cannot make {n} distinct colours from the UCL pool ({len(pool)} available)')
    if n > 1:
        pairs = itertools.combinations(range(len(chosen)), 2)
        closest = min(np.linalg.norm(lab[chosen[i]] - lab[chosen[j]]) for i, j in pairs)
        if closest < min_distance:
            print(f'WARNING: {n} colours requested; the closest two are only {closest:.1f} apart in CIELAB')
    return colours


def categorical_palette(categories: Iterable[str]) -> Dict[str, str]:
    """Map each category to a different UCL-derived colour (same input order gives the same colours)."""
    categories = list(categories)
    return dict(zip(categories, ucl_palette(len(categories))))


def age_palette(categories: Iterable[str]) -> Dict[str, str]:
    """Ordered blue -> purple -> pink -> orange palette for developmental age (categories in age order)."""
    categories = list(categories)
    cmap = mpl.colors.LinearSegmentedColormap.from_list('ucl_age', [HERITAGEBLU, BLU, BRIGHTPUR, FUS, ORA])
    points = np.linspace(0, 1, max(len(categories), 2))
    return {c: mpl.colors.to_hex(cmap(p)) for c, p in zip(categories, points)}


def point_size(n_cells: int, width_in: float) -> float:
    """Scatter point area (pt^2) for a panel `width_in` inches wide. Scanpy's default (120000 / n_cells) is meant for a
    4 inch panel, so it is scaled by panel area; floored so points stay visible for very large datasets."""
    return max(0.3, 120000 / max(n_cells, 1) * (width_in / 4) ** 2)


def usage_cmap():
    """Grey to UCL bright purple, for GEP usage and gene expression."""
    return sb.color_palette(f'blend:gainsboro,{BRIGHTPUR}', as_cmap=True)


def set_publication_style(dpi: int = 300) -> None:
    """Journal-style defaults: small sans-serif text, thin lines, editable text in PDF/SVG, rasterised scatter points."""
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
        'font.size': 7,
        'axes.titlesize': 8,
        'axes.labelsize': 7,
        'xtick.labelsize': 6,
        'ytick.labelsize': 6,
        'legend.fontsize': 6,
        'axes.linewidth': 0.6,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': False,
        'pdf.fonttype': 42,      # real (editable) text in PDFs
        'svg.fonttype': 'none',  # real text in SVGs
        'savefig.dpi': dpi,
        'savefig.transparent': False,
        'savefig.bbox': 'tight',
        'figure.dpi': 100,
    })
    # scanpy scatter points are drawn as an image inside vector files, so PDFs stay small
    sc.set_figure_params(vector_friendly=True, dpi_save=dpi, frameon=False, fontsize=7)
    mpl.rcParams['font.size'] = 7
