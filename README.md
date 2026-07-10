# drg-sn-atlas

Core single-cell RNA-seq analysis code for the integrated **mouse dorsal root ganglion (DRG) and sciatic nerve cell atlas**.

> Meriau P, Thomsen MB, Avraham O, Trauterman B, Feng R, Ewan EE, Cavalli V (2026).
> *Single-cell transcriptomic atlas of glial cells in adult mouse dorsal root ganglia identifies multipotent progenitors.* (bioRxiv — DOI pending)

This repository contains the core pipeline used to build and annotate the atlas: integration of >200,000 DRG and sciatic nerve transcriptomes across multiple published and lab-generated datasets, hierarchical clustering, cell-type annotation, differential expression and GO enrichment, gene co-expression network analysis (WGCNA), cell–cell communication, pseudotime, and QC/export.

> **Status:** core analysis pipeline. Supplementary analyses (reference mapping / label transfer, RNA velocity, injury time-course dynamics) and figure-panel scripts will be added prior to publication.

## Repository structure

```
code/
├── 00_setup.R … 20_main_qc_export.R    # main pipeline, run in numeric order
│      00–05   setup, load & integrate, initial clustering, marker scoring, naming
│      06–07   per-compartment clustering (glia/neurons/immune/fibro/mural/VEC) + label transfer to the full object
│      08–11   UMAPs & module scores, QC plots, DGE + GO, correlations, UpSet plots
│      12–19   glial-compartment analyses: correlations, DGE/GO, WGCNA, CellChat, dotplots, pseudotime, heatmaps
│      20      QC metrics export
├── utils/                              # sourced by every script: libraries, helper functions, palettes/variables, gene lists
└── imports/                            # per-dataset loaders for each integrated published dataset
```

Each script sets its working directory to its own location and uses paths relative to a project root (e.g. `../r_objects/`, `../plots/`). Large Seurat objects and generated outputs are not tracked here — see **Data availability**.

## Software environment

- R 4.3.3, Seurat v5.1.0 (harmony integration)
- Downstream: Monocle3, velocyto.R, speckle (propeller), WGCNA, CellChat, clusterProfiler, LISI
- The full project pins dependencies with `renv`; `code/utils/libraries.R` lists the load-time packages.

## Data availability

- **Newly generated (GEO):** integrated atlas samples — [GSE317728](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE317728); Fmr1-KO DRG — GSE337371; sciatic nerve crush time-course — GSE337372.
- **Processed atlas object (Seurat):** Zenodo — DOI pending.
- **Previously published datasets** reanalysed here are available under their original GEO accessions (listed in the manuscript; loaders in `code/imports/`).

## Citation

If you use this code, please cite the manuscript above (DOI pending) and this repository.

## License

Released under the [MIT License](LICENSE).

## Contact

Cavalli Lab, Washington University in St. Louis.
