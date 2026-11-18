# Spyo-GENOMiCUS-399

Analysis code for **"Pangenomic decomposition reveals lineage-specific immune-evasion
and virulence architectures in group A *Streptococcus*"**
(Chauhan, Monk, Palsson & Nizet).

Non-negative matrix factorization of the *Streptococcus pyogenes* accessory genome across
399 complete genomes, decomposing it into 11 **phylons** — co-inherited accessory gene
modules corresponding to lineage-specific repertoires and independently circulating
prophages.

## What this repository is for

The code here reproduces every figure and table in the paper from the deposited matrices.
Issues and pull requests are welcome: if a step does not reproduce, or you want to apply
the workflow to another species, please open an issue.

## Data

Genome assemblies are public and are **not** redistributed here. `data/accessions.tsv`
lists the 399 BV-BRC and NCBI RefSeq accessions.

All derived data — the pangenome presence/absence matrix **P**, the core/accessory/rare
partitions, the **L** and **A** matrices, and Supplementary Tables S1–S12 — are archived
on Zenodo:

> **DOI: [10.5281/zenodo.20403112](https://doi.org/10.5281/zenodo.20403112)**

That is the *concept* DOI and always resolves to the newest version. Download it with
`scripts/zenodo_fetch.py` and unpack into `data/`.

## Pipeline

Run in order. Steps 1–3 rebuild the pangenome from raw assemblies and take days of compute
on a large server; **steps 4 onward run in minutes from the Zenodo matrices**, which is the
path most readers will want.

| Step | Notebook / script | Produces |
|---|---|---|
| 1a | `scripts/1a_p1_fetch_ncbi_metadata.py`, `1a_p2_filter_ncbi_metadata.py` | NCBI metadata |
| 1a | `notebooks/1a_filter_genomes_for_download_BV-BRC.ipynb` | candidate genome list |
| 1b | `notebooks/1b_download_genomes_BV-BRC.ipynb` | assemblies |
| 2a | `notebooks/2a_clean_metadata.ipynb` | unified metadata table |
| 2b | `notebooks/2b_mash_clustering_and_filtration.ipynb` | Mash QC, species-boundary filter |
| 2c | `notebooks/2c_build_cds_pangenome.ipynb` | CD-HIT and PanTA pangenomes |
| 2d | `notebooks/2d_enrich_metadata.ipynb` | MLST, *emm* types, source categories |
| 3a | `notebooks/3a_extract_CAR_genomes.ipynb` | core / accessory / rare partitions |
| 4a | `notebooks/4a_NMF_decomposition.ipynb` | **L** and **A** matrices, rank calibration |
| 5a | `notebooks/5a_phylon_characterization.ipynb` | phylon naming, Figs 1–6, Table 1 |
| 6a | `notebooks/6a_poppunk_concordance.ipynb` | PopPUNK concordance, Tables S10–S12, Fig S8 |

Helper scripts in `scripts/` are called by the notebooks or run standalone:

- `extract_bakta_annotations.py` — AMR, VFDB, IS and pilin-locus gene content from Bakta JSON
- `make_figure_s8.py` — UpSet plot and zero-phylon frequency panel
- `zenodo_fetch.py` — download and checksum-verify the Zenodo deposit

## Installation

```bash
conda env create -f environment.yml
conda activate spyo-genomicus
```

External tools not installable via conda — Bakta (v1.12.0, database v6.0) and PanTA — are
needed only for steps 1–3. See the Methods section of the paper for their invocation.

## Citing

Please cite the paper and, if you use the deposited matrices, the Zenodo DOI. See
`CITATION.cff`.

## License

MIT — see `LICENSE`.
