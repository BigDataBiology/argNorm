# argNorm

Fast ARG normalization by mapping to the ARO ontology.

## Gene name matching

|Software|Database|Gene name column|Matches with (n)th field of fasta info|Note|
|:-------|:-------|:---------------|:-------------------------------------|:---|
|Abricate|ncbi|`gene_symbol`|5 or 6|`sep='\|'` |
|Abricate|argannot|`gene_symbol` or `gene_name`|2|`sep='~~~'`|
|Abricate|card|`gene_symbol`|4 (in faa, not fna) |`sep='\|'`, some output gene names are prefixed by taxonomy names (e.g. Pseudomonas_aeruginosa_emrE)|
|Abricate|megares|`gene_name`|Exact match|Drop first field and `replace(':', '\|')`|
|Abricate|resfinder|`gene_symbol`|1|`sep='_'` but gene name contains `sep`|
|deeparg|deeparg|`gene_name`|Exact match||
|ARGs-OAP|sarg|||[Known issue](https://github.com/AdeBC/quick_amr_db_harmonisation/issues/2)|