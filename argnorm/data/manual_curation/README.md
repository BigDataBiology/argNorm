# Source

- These manual cuation tables are generated in `db_harmonisation/manual_curation` after running `db_harmonisation/crude_db_harmonisation.py`.
- The script will generate most genes that should be manually curated. These are manually processed, curated, and added to this folder.
- Resfinder gene clusters aren't flagged automatically and are manually collated in resfinder mappings.
- Genes are manually curated by searching them up online on CARD, NCBI GenBank, UniProt, etc.
- Description columns are filled out manually.
- Please refer to `db_harmonisation/README.md` for more information on the automatic flagging of loose RGI hits for manual curation.

# Resfinder Notes

## Gene Clusters

- Resfinder has gene clusters (nucleotide sequence with multiple CDSs present) which can't be passed through RGI using 'contig' mode.
- Gene clusters were identified and were manually assigned ARO numbers.
- 40 gene clusters present.

## Reverse Complement
1) blaBIM-1_1_CP016446
2) blaSPG-1_1_KP109680
3) grdA_1_QJX10702
4) tet(43)_1_GQ244501
5) aac(3)-Xa_1_AB028210
6) blaBKC-1_1_KP689347
7) mph(A)_1_D16251
8) qepA1_1_AB263754
9) aac(3)-I_1_AJ877225

- 9 genes in reverse complement form also present.
- RC genes were manually curated

# NCBI Notes

NCBI also has stress response genes and virulence genes not covered by CARD. This is the reason for the missing mappings in NCBI curation. 