# Data

Genome assemblies are not redistributed here; `accessions.tsv` lists all 399 BV-BRC and
NCBI RefSeq accessions used in the study.

All derived matrices and supplementary tables are archived on Zenodo under the concept DOI
**10.5281/zenodo.20403112**, which always resolves to the newest version.

    python3 ../scripts/zenodo_fetch.py 20403112 -o .

Expected contents after unpacking:

| File | Description |
|---|---|
| `P_matrix_panta.csv` | full pangenome gene presence/absence (8,864 x 399) |
| `df_acc_panta_cdhit.pickle` | accessory-genome matrix (1,167 x 399) |
| `df_rare_panta_cdhit.pickle` | rare-genome matrix (6,405 x 399) |
| `L_binarized_panta_cdhit.csv` | binarized L matrix (1,167 genes x 11 phylons) |
| `A_binarized_panta_cdhit.csv` | binarized A matrix (11 phylons x 399 genomes) |
| `Spyo_Supplementary_Tables.xlsx` | Supplementary Tables S1-S12 |
