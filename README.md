# Codon-usage clustering pipeline

This repository contains a small Python pipeline to (i) cluster genes based on per-gene codon/AA usage profiles, (ii) export results to Excel, and (iii) generate 2D embedding overlay plots for curated gene sets.

## Quick start

1. Install dependencies (example with pip):

```bash
pip install -r requirements.txt
```

2. Run the main pipeline:

```bash
python Clustering_Pipeline.py
```

You will be prompted to select:
- a codon-usage Excel workbook (per-gene counts/usage)
- a GeneIDs/annotation Excel workbook (used for gene symbols/descriptions; optional but recommended)
- a curated cluster file (TXT/CSV/TSV/Excel with one column per cluster, containing locus tags)
- a CDS FASTA file (for identifier canonicalization and sequence-derived metrics)

The pipeline writes:
- a main Excel workbook (coordinates, membership tables, metrics, and usage tables)
- a per-cluster gene list workbook (one sheet per cluster)
- optional cluster-wise codon usage workbooks (ACU/RCU/RCU_devZ), depending on settings

3. Generate plots (optional)

If `PIPELINE['auto_run_plotting_pipeline'] = True`, the pipeline will call `Plotting_Pipeline.py` automatically.
You can also run the plotting script directly:

```bash
python Plotting_Pipeline.py
```

## Notes

- UMAP is optional. If you select `dimred_method='umap'`, install `umap-learn`.
- K-medoids is optional. If `sklearn-extra` is not available, the code falls back to K-means.
- Custom colormaps are supported via an optional Excel workbook, but are disabled by default in the clean configuration.

## Repository layout

- `Clustering_Pipeline.py` — main analysis/export pipeline
- `Plotting_Pipeline.py` — cluster overlay figures from the output workbook
- `codonpipe/` — reusable modules (I/O, metrics, plotting, statistics)

## Colormaps

This repository includes a built-in MATLAB-style **parula** colormap (256 steps).
You can use it by setting any colormap option to `parula` or `parula_r`.
