# Data

Genome assemblies are public and are not redistributed here.

| File | Contents |
|---|---|
| `accessions.tsv` | All 399 genomes with BV-BRC identifiers and NCBI assembly, RefSeq, GenBank, BioProject and BioSample accessions |
| `car_metadata_Sypo.csv` | Curated per-genome metadata (399 x 137): MLST, *emm* type and cluster, isolation source, country, collection date, CheckM metrics |

Everything else the scripts read — the pangenome matrix **P**, the
core/accessory/rare partitions, the binarized **L** and **A** matrices, the PanTA
gene annotations and the Supplementary Tables — is archived on Zenodo under the
concept DOI **10.5281/zenodo.20403112**, which always resolves to the newest
version:

    python scripts/zenodo_fetch.py 20403112 -o data

Expected after unpacking:

| File | Description |
|---|---|
| `L_bin_panta_cdhit_Spyo.csv` | binarized L matrix (1,167 genes x 11 phylons) |
| `A_bin_panta_cdhit_Spyo.csv` | binarized A matrix (11 phylons x 399 genomes) |
| `df_acc_panta_cdhit.pickle` | accessory P matrix (1,167 x 399) |
| `df_genes_panta_cdhit_final.pickle` | full pangenome P matrix (8,864 x 399) |
| `gene_annot_panta_Spyo.csv` | PanTA gene-family annotations |
| `sda1_hits.tsv` | tBLASTn Sda1 ortholog hits |
| `Spyo_Supplementary_Tables.xlsx` | Supplementary Tables S1-S12 |

Downloaded matrices are gitignored; `accessions.tsv` and `car_metadata_Sypo.csv`
are tracked.
