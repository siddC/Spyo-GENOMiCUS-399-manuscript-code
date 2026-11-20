# Spyo-GENOMiCUS-399

Analysis code for **"Pangenomic decomposition reveals lineage-specific immune-evasion and
virulence architectures in group A *Streptococcus*"** — Chauhan, Monk, Palsson & Nizet.

Non-negative matrix factorization of the *Streptococcus pyogenes* accessory genome across
399 complete genomes decomposes it into 11 **phylons**: co-inherited accessory gene modules
corresponding to lineage-specific repertoires and to independently circulating prophages
("mobilons").

Issues and pull requests are welcome. If a step does not reproduce, or you want to apply the
workflow to another species, please open an issue.

## Data

Genome assemblies are public and are **not** redistributed here. `data/accessions.tsv` lists
all 399 genomes with their BV-BRC identifiers and NCBI assembly, BioProject and BioSample
accessions; `data/car_metadata_Spyo.csv` is the full curated metadata table used throughout.

Derived data — the pangenome matrix **P**, the core/accessory/rare partitions, the **L** and
**A** matrices, and Supplementary Tables S1–S12 — are archived on Zenodo:

> **[10.5281/zenodo.20403112](https://doi.org/10.5281/zenodo.20403112)**

That is the *concept* DOI and always resolves to the newest version.

```bash
python scripts/zenodo_fetch.py 20403112 -o data
```

## Installation

```bash
conda env create -f environment.yml
conda activate spyo-genomicus
```

Bakta (v1.12.0, database v6.0) and PanTA are not installable via conda and are needed only
for the pangenome construction steps; see the paper's Methods for their invocation.

## Running the code

All scripts read from `data/` and write to `results/` when run **from the repository root**.
Either directory can be overridden:

```bash
python scripts/table1_analysis.py                      # defaults
SPYO_DATA=/mnt/craig/… python scripts/table1_analysis.py
python scripts/phylon_and_global_concordance.py /path/to/data /path/to/results
```

Steps 1–3 rebuild the pangenome from raw assemblies and take days of compute on a large
server. **Steps 4 and 5 run in minutes from the Zenodo matrices**, which is the path most
readers will want.

## Notebooks

Organised by OSEMN stage under `notebooks/S_pyogenes/`.

| Stage | Notebook | Produces |
|---|---|---|
| **1 — Obtain** | `1_obtain_data/1a_filter_genomes_for_download_BV-BRC.ipynb` | candidate genome list from BV-BRC |
| | `1_obtain_data/1b_download_genomes_BV-BRC.ipynb` | assemblies |
| **2 — Scrub** | `2_scrub_data/2a_clean_metadata.ipynb` | unified metadata table |
| | `2_scrub_data/2b_mash_clustering_and_filtration_Spyogenes.ipynb` | Mash QC, species-boundary filter |
| | `2_scrub_data/2c_build_cds_pangenome.ipynb` | CD-HIT and PanTA pangenomes |
| | `2_scrub_data/2d_enrich_metadata.ipynb` | MLST, *emm* types, source categories |
| **3 — Explore** | `3_explore_data/3a_extract_CAR_genomes.ipynb` | core / accessory / rare partitions |
| | `3_explore_data/final_gene_list.ipynb` | final gene-family list |
| **4 — Model** | `4_model_data/4a_NMF_decomposition.ipynb` | rank calibration, **L** and **A** matrices |
| **5 — Interpret** | `5_interpret_data/5a_phylon_characterization.ipynb` | phylon naming, Figs 1–6, Table 1 |
| | `5_interpret_data/5b_poppunk_concordance_and_amr.ipynb` | PopPUNK concordance, AMR enrichment, Tables S10–S12, Fig. S8 |
| | `5_interpret_data/sda1_extraction.ipynb`, `Sda1_prophage_extraction.ipynb` | Sda1 prophage locus extraction |
| | `5_interpret_data/spec_mf2_extraction.ipynb` | *speC*/*mf2* cassette extraction |
| | `5_interpret_data/gene_pa_core_genome_analysis.ipynb` | core-genome occupancy analysis |
| | `5_interpret_data/A_norm_best_genome.ipynb` | representative-genome selection by A-matrix affinity |
| | `5_interpret_data/Untitled.ipynb` | miscellaneous supporting calculations |

Helper scripts alongside the notebooks drive the external tools for stages 1–2:
`download_bvbrc_metadata.py`, `run_bakta.py`, `run_mash.sh`, `run_checkm.py`,
`merge_genome_metadata.py`, and under `1_obtain_data/ncbi/`:
`1a_p1_fetch_ncbi_metadata.py`, `1a_p2_filter_ncbi_metadata.py` and
`download_ncbi_genomes.py`.

## Scripts

Standalone reproducibility scripts in `scripts/`, each regenerating a specific published
result.

| Script | Reproduces |
|---|---|
| `phylon_and_global_concordance.py` | Per-phylon and global concordance vs *emm* and MLST (Tables S5, S6) |
| `table1_analysis.py` | Table 1 — per-phylon virulence and ecological characterization |
| `sda1_enrichment_analysis.py` | Table S7 — Sda1 prophage core enrichment (Fisher, BH-corrected) |
| `sda1_results_verification.py` | Verification of every Sda1 claim in the Results |
| `nested_sda1_prophage_analysis.py` | Table S8 — nested prophage in strain 1314.4589 |
| `speC_mf2_synteny_table.py` | Table S9 — *speC*/*mf2* synteny annotations |
| `regulon_emmswitching_verification.py` | Verification of the *emm*-pattern regulon and *emm*-switching results |
| `extract_bakta_annotations.py` | AMR, VFDB, IS and pilin-locus gene content from Bakta JSON |
| `fig_sda1_prophage_synteny.py` | Figure 5d — Sda1 prophage cross-phylon synteny |
| `fig_speC_mf2_synteny.py` | Figure 6 — *speC*/*mf2* cassette synteny |
| `figS_double_sda1_prophage.py` | Figure S6 — double Sda1 prophage |
| `make_figure_s8.py` | Figure S8 — UpSet plot and zero-phylon frequency distribution |
| `zenodo_fetch.py` | Download and checksum-verify the Zenodo deposit |

### PopPUNK

Population-structure assignment is run outside the notebook (PopPUNK v2.7.8, reference
database `Streptococcus_pyogenes_v1_refs`; Davies et al., *Nat Genet* 2019). The exact
commands, including the working-directory caveat that avoids PopPUNK's plotting bug on
`--output` paths containing a separator, are documented in
`5_interpret_data/5b_poppunk_concordance_and_amr.ipynb`.

## Citing

Please cite the paper and, if you use the deposited matrices, the Zenodo DOI. See
`CITATION.cff`.

## License

MIT — see `LICENSE`.
