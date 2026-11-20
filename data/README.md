# Data

Genome assemblies are public and are not redistributed here.

## Tracked in this repository

| File | Contents |
|---|---|
| `accessions.tsv` | All 399 genomes with BV-BRC identifiers and NCBI assembly, RefSeq, GenBank, BioProject and BioSample accessions |
| `car_metadata_Spyo.csv` | Curated per-genome metadata (399 x 137): MLST, *emm* type and cluster, isolation source, country, collection date, CheckM metrics |
| `GAS_assigned_clusters.csv` | PopPUNK genomic-phylogroup assignments against the reference database (399 genomes in 104 clusters, four of which merge two reference phylogroups) |
| `GAS_refined_clusters.csv` | PopPUNK *de novo* cluster assignments from the refined fit (399 genomes in 225 clusters) |

The two PopPUNK tables are small and are kept here so that
`5_interpret_data/5b_poppunk_concordance_and_amr.ipynb` runs without a download step.

## Archived on Zenodo

Everything else the scripts read is archived under the *concept* DOI
**10.5281/zenodo.20403112**, which always resolves to the newest version:

    python scripts/zenodo_fetch.py 20403112 -o data

Expected after unpacking:

| File | Description |
|---|---|
| `L_bin_panta_cdhit_Spyo.csv` | binarized L matrix (1,167 accessory genes x 11 phylons) |
| `A_bin_panta_cdhit_Spyo.csv` | binarized A matrix (11 phylons x 399 genomes) |
| `df_core_panta_cdhit.pickle` | core P matrix (1,292 x 399) |
| `df_acc_panta_cdhit.pickle` | accessory P matrix (1,167 x 399) |
| `df_rare_panta_cdhit.pickle` | rare P matrix (6,405 x 399) |
| `gene_annot_panta_Spyo.csv` | curated PanTA gene-family annotations |
| `Spyo_bakta_annot_expert_long.tsv` | Bakta expert-system calls, long form (AMRFinderPlus, VFDB, IS, BlastRules) |
| `Spyo_bakta_annot_fct_long.tsv` | Bakta FCT / pilin-locus gene calls, long form |
| `sda1_hits.tsv` | tBLASTn Sda1 ortholog hits |
| `Spyo_Supplementary_Tables.xlsx` | Supplementary Tables S1-S12 |

The record also holds the raw PanTA and CD-HIT outputs (`gene_presence_absence.csv`,
`gene_annotation.csv`, the representative-sequence FASTAs, the Mash and BLAST
intermediates) and the stage 1-4 notebooks. Those are not needed to reproduce the
published tables and figures from stage 5.

Three files appear twice under different names. `L_binarized_panta_cdhit.csv`,
`A_binarized_panta_cdhit.csv` and `Spyo_GENOMiCUS_399_Suppl_Tables.xlsx` are
byte-identical copies retained under the names used in earlier versions of the
record. Use the names in the table above.

### Reconstructing the full pangenome matrix

`scripts/sda1_enrichment_analysis.py` and `scripts/sda1_results_verification.py`
expect `df_genes_panta_cdhit_final.pickle`, the full 8,864 x 399 matrix. It is not
deposited separately because the three partitions above reconstruct it exactly —
they share identical columns and have no overlapping rows:

    import pandas as pd

    df = pd.concat([
        pd.read_pickle("data/df_core_panta_cdhit.pickle"),   # 1,292 x 399
        pd.read_pickle("data/df_acc_panta_cdhit.pickle"),    # 1,167 x 399
        pd.read_pickle("data/df_rare_panta_cdhit.pickle"),   # 6,405 x 399
    ])
    assert df.shape == (8864, 399)
    df.to_pickle("data/df_genes_panta_cdhit_final.pickle")

## What is and is not under version control

`data/*.csv`, `data/*.pickle` and `data/*.xlsx` are gitignored so that a Zenodo
fetch does not dirty the working tree. The four files in the first table above are
explicit exceptions in `.gitignore` and are tracked.
