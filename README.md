# drg-sn-atlas

Analysis code for the integrated mouse dorsal root ganglion (DRG) and sciatic
nerve single-cell RNA-seq atlas (Meriau, Thomsen et al., 2026). The atlas
integrates >200,000 cells across multiple published and lab-generated datasets
and annotates them at three hierarchical levels (Class, Major, Minor).

## Repository structure

```
code/
  00_setup.R … 20_main_qc_export.R   main pipeline, run in numeric order
    00–05   setup, loading and integration, initial clustering, marker scoring, naming
    06–07   per-compartment clustering (glia, neurons, immune, fibroblast, mural, vascular)
            and label transfer back to the full object
    08–11   UMAPs and module scores, QC plots, DGE and GO, correlations, UpSet plots
    12–19   glial-compartment analyses: correlations, DGE/GO, WGCNA, CellChat, dotplots,
            pseudotime, heatmaps
    20      QC metrics export
  utils/    shared functions, variables, palettes and gene lists sourced by the scripts
  imports/  per-dataset loaders used during integration
```

Scripts are run from `code/` and read and write sibling data directories
(`../r_objects/`, `../plots/`, `../markers/`) that are not tracked here.

## Software environment

R 4.3.3 with Seurat 5.1.0 (Harmony integration). Downstream analyses use Monocle3,
velocyto.R, speckle, WGCNA, CellChat, clusterProfiler and LISI. Dependencies are
pinned with renv in the full project; the load-time package list is in
`code/utils/libraries.R`.

## Data

- Sequencing data (GEO): GSE317728, GSE337371, GSE337372
- Processed atlas object (Seurat): https://doi.org/10.5281/zenodo.21299570
- Previously published datasets reanalysed here are available under their original
  GEO accessions (loaders in `code/imports/`)

## Citation

Meriau P, Thomsen MB, Avraham O, Trauterman B, Feng R, Ewan EE, Cavalli V (2026).
Single-cell transcriptomic atlas of glial cells in adult mouse dorsal root ganglia
identifies multipotent progenitors.

## License

MIT (see LICENSE).
