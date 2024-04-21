# Gene Clusters

- Resfinder has gene clusters (nucleotide sequence with multiple CDSs present) which can't be passed through RGI using 'contig' mode.
- Gene clusters were identified and were manually assigned ARO numbers.
- A seperate file with manual curation for gene clusters and RCs was created, and their AROs were updated after concatenating RGI results and genes not in RGI results.
- 40 gene clusters present.

# Reverse Complement
1) blaBIM-1_1_CP016446
2) blaSPG-1_1_KP109680
3) grdA_1_QJX10702
4) tet(43)_1_GQ244501


- 4 genes in reverse complement form also present.
- blaBIM-1_1_CP016446 and mph(D)_1_AB048591 were not found in CARD and were given parent ARO mappings.
- RGI correctly assigned ARO numbers to other two.

# erm(X)_1_M36726
- Had one to many ARO mapping due to include loose. Included correct ARO number for it.