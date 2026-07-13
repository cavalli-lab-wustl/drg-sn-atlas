# drg-sn-atlas

Analysis code for the integrated mouse dorsal root ganglion and sciatic nerve
single-cell RNA-seq atlas (Meriau, Thomsen et al., 2026).

## Layout

Pipeline scripts are in `code/`, run in numeric order (`00`–`20`): loading and
integration, per-compartment clustering, label transfer, marker/DGE and GO
analysis, correlations, WGCNA, CellChat, pseudotime, heatmaps, and QC export.
Shared functions, variables and gene lists are in `code/utils/`; the per-dataset
loaders used during integration are in `code/imports/`.

Scripts are run from `code/` and read and write sibling data directories
(`../r_objects/`, `../plots/`, `../markers/`) that are not tracked here.

## Environment

R 4.3.3, Seurat 5.1.0. Full package list in `code/utils/libraries.R`.

## Data

- GEO: GSE317728, GSE337371, GSE337372
- Processed atlas object (Seurat): https://doi.org/10.5281/zenodo.21299570

## License

MIT (see LICENSE).
