# argNorm

Fast ARG normalization by mapping to the ARO ontology.

## Gene name matching

|Software|Database|Gene name column|Matches with (n)th field of fasta info|Note|
|:-------|:-------|:---------------|:-------------------------------------|:---|
|Abricate|ncbi|`gene_symbol`|5 or 6|`sep='|'`|
|Abricate|argannot|`gene_symbol` or `gene_name`|||
|Abricate|card|`gene_symbol`|4 (*faa not fna) |`sep='|'`, some output gene names are prefixed by taxonomy names (e.g. Pseudomonas_aeruginosa_emrE)|
|Abricate|megares|`gene_symbol`|||
|Abricate|resfinder|`gene_symbol`|1|`sep='_'` but gene name contains `sep`|
|deeparg|deeparg||||
|ARGs-OAP|sarg||||