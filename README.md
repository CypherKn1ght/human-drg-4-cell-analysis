# DRG neuron immunofluorescence quantification

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

### **Soma segmentation** 
An Otsu threshold picks out bright cell cores. These are used as seeds only. Touching seeds are separated with a distance-transform watershed.

Each seed then expands into a second, lower threshold (Li), which picks up the dimmer parts of a soma. Dim debris with no seed inside it stays excluded.

The expansion is confined to a region that has been morphologically opened, which removes thin structures. Masks therefore stop at the soma edge and do not follow axons.

Masks are dropped if they touch the frame edge, fall outside the size range, or fall below the solidity cutoff.

### **Neurite segmentation** 
A Sato tubeness filter scores each pixel on how tube-like it is. Thresholding on brightness alone does not work here, because the cutoff needed to catch faint neurites also catches background. Somas are dilated and subtracted, leaving the processes.

Neurites are pooled into one measurement per field. They are not assigned to individual cells, because processes cross and overlap and there is no way to tell which neuron a segment belongs to in a single 2D image. Adding them to cell masks would also pull each cell's mean toward the dimmer neurite value.

### **Measurement**
Background is measured per image and per channel from pixels far from any mask, then subtracted.

Masks are eroded a few pixels before measuring, so blurred edge pixels are excluded.

The script reports mean intensity. Integrated density scales with cell size, and DRG soma diameter varies too much for that to be useful.

### **Correlations** 
Correlations are calculated separately for each plate. Plates often differ in overall intensity, and combining them can hide a relationship that holds within every plate.

A combined value is also reported. It uses z-scores computed within each plate, which removes the offsets.

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
