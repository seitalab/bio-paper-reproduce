# Paper reproduction

Notebooks and notes for reproducing analyses and figures from research papers. This repository contains data download scripts, analysis notebooks, figures, and writeups.

## Setup

From the project root, create the environment with Conda (or replace `conda` with `mamba` in the creation command):

```bash
conda env create --prefix ./.venv --file environment.yml
conda activate "$PWD/.venv"
jupyter lab
```

Open a notebook in `notebooks/` and select the environment's Python kernel. Downloaded inputs are stored in `data/`; figures and writeups are in `figures/` and `docs/`.

## Papers

### A sustained type I IFN-neutrophil-IL-18 axis drives pathology during mucosal viral infection

- **Notebook:** [ifn-neutrophil-reproduce.ipynb](notebooks/ifn-neutrophil-reproduce.ipynb)
- **Writeup:** [Reproduction notes](docs/type-I-IFN-neutrophil-IL-18-reproduce-writeup.md)

Download the data into `data/ifn-neutrophil/` from the project root:

```bash
python scripts/ifn-neutrophil-dataset-download.py
```

#### Acknowledgements

- **Paper:** Lebratti et al. (2021), [A sustained type I IFN-neutrophil-IL-18 axis drives pathology during mucosal viral infection](https://doi.org/10.7554/eLife.65762), *eLife*.
- **Dataset:** The authors' single-cell RNA-seq data, available through GEO as [GSE161336](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161336).
- **Gene set:** [MSigDB Hallmark IFN-alpha response](https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_INTERFERON_ALPHA_RESPONSE).
