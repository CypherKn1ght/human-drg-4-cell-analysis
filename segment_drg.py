"""
DRG neuron immunofluorescence quantification.

Segments neuronal cell bodies on the FITC channel, then measures background-
subtracted Cy5 and Texas Red intensity within those regions. Neurites are
segmented separately and measured as their own compartment.

Pipeline
    1. Load channels. Handles 16-bit TIFF and pseudocolored 8-bit JPG.
    2. Seeded region growing on FITC:
         - Otsu threshold finds bright soma cores (seeds only)
         - distance-transform watershed splits touching cores
         - each core grows outward into a permissive Li threshold, so the
           dim rim of a soma is included while isolated dim debris is not
         - growth is bounded by a morphological opening, so masks stop at
           the soma edge instead of running down an axon
    3. Neurites via Sato tubeness ridge filter, with somas subtracted.
    4. Per-image background from cell-free regions, subtracted per channel.
    5. Masks eroded a few px before measurement to avoid blurry edge pixels.

Statistics note
    Correlations are computed WITHIN each plate. Pooling raw values across
    plates that sit at different offsets destroys a real within-plate
    relationship (Simpson's paradox). A pooled value is reported on
    within-plate z-scores instead.

Usage
    python segment_drg.py "path/to/image/folder"

Expects files named *_fitc.*, *_cy5.*, *_txred.* sharing a common prefix.
Writes an "analysis" subfolder containing per-cell and per-plate CSVs,
QC overlays, and summary figures.
"""

import sys
import os
import glob
import hashlib
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from skimage import io, filters, morphology, measure, segmentation, exposure
from scipy import ndimage as ndi

# --- tuning ---------------------------------------------------------------
MIN_AREA = 1500           # px, smallest accepted soma
MAX_AREA = 60000          # px, largest accepted soma
MIN_SOLIDITY = 0.75       # 1.0 = perfectly convex; filters stringy fragments
SEED_SCALE = 1.00         # x Otsu.     lower -> more dim cells seeded
GROW_SCALE = 1.00         # x Li.       lower -> somas grow fatter
NEURITE_SCALE = 0.40      # x triangle. lower -> fainter processes included
ERODE_PX = 2              # shrink mask before measuring
RIDGE_SIGMAS = (1, 2, 4)  # px half-widths of processes to detect
# --------------------------------------------------------------------------


def _drop_small(mask, min_px):
    """Remove connected components smaller than min_px. Version-proof."""
    lab, n = ndi.label(mask)
    if n == 0:
        return mask
    sizes = np.bincount(lab.ravel())
    sizes[0] = 0
    return np.isin(lab, np.flatnonzero(sizes >= min_px))


def load_gray(path):
    """Return a 2D float array. Handles 16-bit TIFF and pseudocolored JPG."""
    im = io.imread(path)
    if im.ndim == 2:
        return im.astype(float)
    return im.astype(float).max(axis=2)


def segment(fitc):
    """Return (soma_labels, neurite_mask) from the FITC channel."""
    sm = filters.gaussian(fitc, sigma=2, preserve_range=True)

    # seeds: bright cores only
    seed_mask = sm > filters.threshold_otsu(sm) * SEED_SCALE
    seed_mask = ndi.binary_fill_holes(seed_mask)
    seed_mask = morphology.opening(seed_mask, morphology.disk(6))
    seed_mask = _drop_small(seed_mask, MIN_AREA)

    # split touching cores
    dist = filters.gaussian(ndi.distance_transform_edt(seed_mask), sigma=4)
    markers, _ = ndi.label(morphology.local_maxima(dist))
    cores = segmentation.watershed(-dist, markers, mask=seed_mask)

    # grow region: permissive, still soma-shaped
    grow = sm > filters.threshold_li(sm) * GROW_SCALE
    grow = ndi.binary_fill_holes(grow)
    grow = morphology.opening(grow, morphology.disk(10))
    grow = morphology.dilation(grow, morphology.disk(2))
    grow |= seed_mask

    labels = segmentation.watershed(-sm, cores, mask=grow)
    labels = segmentation.clear_border(labels)

    out = np.zeros_like(labels)
    keep = 0
    for r in measure.regionprops(labels):
        if MIN_AREA <= r.area <= MAX_AREA and r.solidity >= MIN_SOLIDITY:
            keep += 1
            out[labels == r.label] = keep

    # neurites: elongated ridges, excluding somas
    ridge = filters.sato(sm, sigmas=RIDGE_SIGMAS, black_ridges=False)
    nmask = ridge > filters.threshold_triangle(ridge) * NEURITE_SCALE
    nmask &= sm > filters.threshold_triangle(sm) * NEURITE_SCALE
    nmask &= ~morphology.dilation(out > 0, morphology.disk(4))
    nmask = _drop_small(nmask, 200)

    return out, nmask


def background(img, occupied):
    """Median intensity well away from any segmented structure."""
    far = ~morphology.dilation(occupied, morphology.disk(25))
    return float(np.median(img[far]))


def save_qc(fitc, somas, neurites, tag, outdir):
    """Write a PNG showing soma outlines (red) and neurite mask (blue)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    g = exposure.rescale_intensity(fitc, in_range=(0, 110)).astype(float)
    g /= max(g.max(), 1)
    rgb = np.dstack([g] * 3)
    rgb[neurites] = [0.15, 0.55, 1.0]
    rgb[segmentation.find_boundaries(somas, mode="thick")] = [1, 0.1, 0.1]

    plt.figure(figsize=(14, 12))
    plt.imshow(rgb)
    plt.axis("off")
    for r in measure.regionprops(somas):
        plt.text(r.centroid[1], r.centroid[0], str(r.label),
                 color="yellow", fontsize=9, ha="center")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"qc_{tag}.png"), dpi=85,
                bbox_inches="tight")
    plt.close()


def measure_field(fitc_path, cy5_path, txred_path, tag, qcdir=None):
    fitc, cy5, txr = (load_gray(p) for p in (fitc_path, cy5_path, txred_path))
    somas, neurites = segment(fitc)
    if somas.max() == 0:
        return pd.DataFrame(), {}
    if qcdir:
        save_qc(fitc, somas, neurites, tag, qcdir)

    occupied = (somas > 0) | neurites
    bg = {"cy5": background(cy5, occupied), "txred": background(txr, occupied)}
    cy5c = np.clip(cy5 - bg["cy5"], 0, None)
    txrc = np.clip(txr - bg["txred"], 0, None)

    eroded = np.zeros_like(somas)
    for lab in range(1, somas.max() + 1):
        eroded[morphology.erosion(somas == lab,
                                  morphology.disk(ERODE_PX))] = lab

    df = pd.DataFrame([
        {"label": r.label, "area_px": r.area,
         "diameter_px": r.equivalent_diameter_area, "solidity": r.solidity,
         "cy5_mean": r.intensity_mean}
        for r in measure.regionprops(eroded, intensity_image=cy5c)
    ]).set_index("label")
    for name, img in (("txred_mean", txrc), ("fitc_mean", fitc)):
        df[name] = pd.Series({
            r.label: r.intensity_mean
            for r in measure.regionprops(eroded, intensity_image=img)
        })
    df["field"] = tag
    df = df.reset_index()

    field = {
        "field": tag,
        "n_cells": len(df),
        "soma_cy5": df.cy5_mean.mean(),
        "soma_txred": df.txred_mean.mean(),
        "neurite_cy5": float(cy5c[neurites].mean()) if neurites.any() else np.nan,
        "neurite_txred": float(txrc[neurites].mean()) if neurites.any() else np.nan,
        "neurite_area_px": int(neurites.sum()),
        "neurite_frac_of_field": float(neurites.mean()),
        "bg_cy5": bg["cy5"],
        "bg_txred": bg["txred"],
        "fitc_saturated_frac": float(np.mean(fitc[somas > 0] >= 254)),
    }
    return df, field


def plate_of(name):
    """Plate identifier: the field name with its trailing _NNNN removed."""
    parts = name.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return name


def fitc_hash(path):
    """Hash of decoded pixel content, used to drop duplicate fields."""
    im = io.imread(path)
    arr = im if im.ndim == 2 else im.max(axis=2)
    return hashlib.md5(np.ascontiguousarray(arr)).hexdigest()


def main(folder):
    outdir = os.path.join(folder, "analysis")
    qcdir = os.path.join(outdir, "qc")
    figdir = os.path.join(outdir, "figures")
    for d in (qcdir, figdir):
        os.makedirs(d, exist_ok=True)

    pats = ["*fitc*.tif", "*fitc*.tiff", "*fitc*.jpg", "*fitc*.png"]
    found = sorted({p for pat in pats
                    for p in glob.glob(os.path.join(folder, pat))})
    if not found:
        sys.exit("No files matching *fitc* in " + folder)

    # pair channels, then collapse pixel-identical fields
    seen, jobs, dropped, incomplete = {}, [], [], []
    for f in found:
        stem = f.replace("fitc", "{}")
        cy5, txr = stem.format("cy5"), stem.format("txred")
        name = os.path.basename(f).rsplit("_fitc", 1)[0]
        if not (os.path.exists(cy5) and os.path.exists(txr)):
            incomplete.append(name)
            continue
        h = fitc_hash(f)
        if h in seen:
            dropped.append((name, seen[h]))
            continue
        seen[h] = name
        jobs.append((f, cy5, txr, name))

    if incomplete:
        print("Excluded (missing a channel): " + ", ".join(incomplete))
    if dropped:
        print("Excluded (pixel-identical duplicate of an earlier field):")
        for d, orig in dropped:
            print(f"   {d}  ==  {orig}")
    print(f"\nAnalysing {len(jobs)} unique field(s)\n")

    cells, fields = [], []
    for f, cy5, txr, name in jobs:
        d, fd = measure_field(f, cy5, txr, name.replace(" ", "_"), qcdir=qcdir)
        if not len(d):
            print(f"{name}: no cells found")
            continue
        plate = plate_of(name)
        d["plate"] = plate
        fd["plate"] = plate
        print(f"{name}: {len(d)} cells")
        cells.append(d)
        fields.append(fd)

    if not cells:
        sys.exit("No cells measured.")

    cells = pd.concat(cells, ignore_index=True)
    pd.DataFrame(fields).to_csv(os.path.join(outdir, "per_field.csv"),
                                index=False)

    summ = cells.groupby("plate").agg(
        n_cells=("cy5_mean", "size"),
        cy5_mean=("cy5_mean", "mean"), cy5_sd=("cy5_mean", "std"),
        txred_mean=("txred_mean", "mean"), txred_sd=("txred_mean", "std"),
        diam_mean=("diameter_px", "mean"), diam_sd=("diameter_px", "std"),
    ).round(2)
    summ["cy5_cv"] = (summ.cy5_sd / summ.cy5_mean).round(3)
    summ["txred_cv"] = (summ.txred_sd / summ.txred_mean).round(3)
    summ.to_csv(os.path.join(outdir, "summary_by_plate.csv"))

    # correlations within plate; see module docstring
    percorr = pd.DataFrame([{
        "plate": p, "n": len(g),
        "cy5_vs_txred": g.cy5_mean.corr(g.txred_mean, method="spearman"),
        "cy5_vs_diam": g.cy5_mean.corr(g.diameter_px, method="spearman"),
        "txred_vs_diam": g.txred_mean.corr(g.diameter_px, method="spearman"),
        "cy5_vs_fitc": g.cy5_mean.corr(g.fitc_mean, method="spearman"),
        "txred_vs_fitc": g.txred_mean.corr(g.fitc_mean, method="spearman"),
    } for p, g in cells.groupby("plate")]).round(3)
    percorr.to_csv(os.path.join(outdir, "correlations_by_plate.csv"),
                   index=False)

    z = lambda x: (x - x.mean()) / x.std()
    cells["cy5_z"] = cells.groupby("plate").cy5_mean.transform(z)
    cells["txred_z"] = cells.groupby("plate").txred_mean.transform(z)
    r_spear = cells.cy5_z.corr(cells.txred_z, method="spearman")
    r_naive = cells.cy5_mean.corr(cells.txred_mean, method="spearman")
    cells.to_csv(os.path.join(outdir, "cells.csv"), index=False)

    print("\n" + summ.to_string())
    print("\nSpearman correlations WITHIN each plate:")
    print(percorr.to_string(index=False))
    print(f"\nCy5 vs TxRed, pooled on within-plate z-scores: rho = {r_spear:.3f}")
    print(f"   (naive pooling of raw values gives {r_naive:.3f}, which is")
    print("    misleading when plates sit at different offsets)")

    hi = percorr.cy5_vs_fitc.abs().max()
    if hi > 0.8:
        print(f"\n** WARNING: Cy5 tracks FITC at rho up to {hi:.2f}. That is too")
        print("   tight for independent stains and suggests bleed-through or a")
        print("   thickness/focus artifact. Run a FITC-only control imaged in")
        print("   the Cy5 channel before interpreting any Cy5 result.")

    for p, g in cells.groupby("plate"):
        cut = g.diameter_px.median() * 1.5
        bad = g[(g.diameter_px > cut) | (g.solidity < 0.86)]
        if len(bad):
            print(f"\n** {p}: {len(bad)} cell(s) look merged or misshapen "
                  f"(labels {sorted(bad.label.tolist())}). Check QC overlay.")

    # figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plates = sorted(cells.plate.unique())
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(plates), 3)))

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for a, col, lab in zip(ax, ["cy5_mean", "txred_mean"], ["Cy5", "TxRed"]):
        for p, c in zip(plates, colors):
            a.hist(cells.loc[cells.plate == p, col], bins=12, alpha=0.55,
                   label=p, color=c)
        a.set_xlabel(f"{lab} intensity per soma (a.u.)")
        a.set_ylabel("cells")
        a.legend(fontsize=8)
        a.set_title(f"{lab} distribution")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "distributions.png"), dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    for p, c in zip(plates, colors):
        sub = cells[cells.plate == p]
        ax[0].scatter(sub.cy5_z, sub.txred_z, s=34, alpha=.75, color=c, label=p)
        ax[1].scatter(sub.diameter_px, sub.cy5_mean, s=34, alpha=.75,
                      color=c, label=p)
    ax[0].set_xlabel("Cy5 (z within plate)")
    ax[0].set_ylabel("TxRed (z within plate)")
    ax[0].set_title(f"Co-expression per cell (within-plate z, rho={r_spear:.2f})")
    ax[1].set_xlabel("soma diameter (px)")
    ax[1].set_ylabel("Cy5 (a.u.)")
    ax[1].set_title("Intensity vs soma size")
    for a in ax:
        a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "scatter.png"), dpi=130)
    plt.close(fig)

    print("\nOUTPUT FOLDER: " + os.path.abspath(outdir))
    print("  cells.csv  per_field.csv  summary_by_plate.csv")
    print("  correlations_by_plate.csv")
    print("  qc/*.png       <- check these first")
    print("  figures/*.png")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
