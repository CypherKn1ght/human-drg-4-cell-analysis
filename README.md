# DRG neuron immunofluorescence quantification

Per-cell intensity quantification for cultured human dorsal root ganglion
neurons imaged in three fluorescence channels. Cell bodies are located on a
pan-neuronal FITC stain, then Cy5 and Texas Red intensities are measured within
those regions. Neurites are segmented separately and measured as their own
compartment.

## Install

```bash
pip install -r requirements.txt
```

Python 3.9+, with scikit-image, numpy, pandas, scipy and matplotlib.

## Run

```bash
python segment_drg.py "path/to/images"
```

Expects one file per channel per field, sharing a prefix:

```
ngf plate 1_0000_fitc.jpg
ngf plate 1_0000_cy5.jpg
ngf plate 1_0000_txred.jpg
```

- Plate identity is the field name with its trailing `_NNNN` removed.
- TIFF, PNG and JPG all work. 16-bit TIFF is preferred; 8-bit JPG is lossy and
  clips bright cells at the sensor ceiling.
- Fields missing a channel are skipped. Fields whose FITC image is
  pixel-identical to an earlier one are dropped.

## Output

Written to `analysis/` inside the image folder.

| File | Contents |
|---|---|
| `cells.csv` | One row per neuron: area, diameter, solidity, background-subtracted Cy5 and Texas Red, raw FITC, within-plate z-scores |
| `summary_by_plate.csv` | n, mean, SD and coefficient of variation per plate |
| `correlations_by_plate.csv` | Within-plate Spearman correlations |
| `per_field.csv` | Field-level values, neurite compartment, background, saturated fraction |
| `qc/` | Segmentation overlays: red = somas, blue = neurites, yellow = cell IDs |
| `figures/` | Intensity distributions and per-cell scatter |

Check the overlays before using the numbers. The `label` column in `cells.csv`
matches the numbers drawn on each cell.

The script prints warnings when a measured channel correlates with FITC above
Spearman rho 0.8, when a mask's size or shape suggests two cells were merged,
and when soma pixels sit at the sensor ceiling.

## How it works

**Soma segmentation** — seeded region growing, not a single threshold. An Otsu
threshold finds bright cell cores, used only as seeds; touching cores are split
by a distance-transform watershed. Each seed then grows outward into a more
permissive Li threshold, so the dim side of a soma is captured while isolated
dim debris, having no seed, is not. A morphological opening bounds the growth
region, so masks stop at the soma edge instead of running down an axon. Objects
touching the frame edge, outside the accepted size range, or below the
convexity floor are discarded.

**Neurite segmentation** — a Sato tubeness filter, which scores each pixel on
how much it resembles an elongated tube rather than a blob or noise. A
threshold alone cannot isolate processes, since any cutoff low enough to catch
faint neurites also catches background fluctuation. Somas are dilated and
subtracted from the result.

Neurites are measured as one field-level compartment, not assigned to
individual cells. Processes from different neurons cross and overlap, so
per-cell attribution is not recoverable from a single 2D field, and merging
them into cell masks would dilute soma signal in proportion to how much neurite
surrounded each cell.

**Measurement** — background is estimated per image and per channel from
regions well away from any segmented structure, then subtracted. Each mask is
eroded a few pixels so blurry edge pixels do not contribute. Mean intensity is
reported rather than integrated density, since soma diameter varies widely and
integrated density would scale with cell size.

**Correlations** — computed within each plate, not across the pooled set.
Plates commonly sit at different intensity offsets, and pooling raw values
across them can obscure a relationship present in every plate individually. A
pooled figure is also reported, computed on within-plate z-scores so the
offsets cancel.

## Caveats

- **Cross-channel intensities are not comparable.** Different fluorophores,
  antibodies and exposure times mean each channel carries its own arbitrary
  units. Compare within a channel across conditions, or use rank correlation
  across channels.
- **Thresholds are recomputed per image.** A dense field gets a different
  cutoff than a sparse one. Acceptable within a session with matched exposure;
  consider a fixed threshold for cross-condition work.
- **Diameters are in pixels.** Convert using the objective and camera pixel
  size.
- **Cells within a field are not independent replicates** for plate-level
  questions.

Segmentation parameters are the constants at the top of `segment_drg.py`.
Change one at a time and re-check the overlays.
