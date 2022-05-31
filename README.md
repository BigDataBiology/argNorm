# argNorm

Fast ARG normalization by mapping to the ARO ontology.

## Gene name matching

|Software|Database|Gene name column|Matches with (n)th field of fasta info|Note|
|:-------|:-------|:---------------|:-------------------------------------|:---|
|[ABRicate](https://github.com/tseemann/abricate)|[ncbi](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt)|`gene_symbol`|5 or 6|`sep='\|'` |
|ABRicate|[argannot](https://github.com/tseemann/abricate/tree/master/db/argannot)|`gene_symbol` or `gene_name`|2|`sep='~~~'`|
|ABRicate|card|`gene_symbol`|4 (in faa, not fna) |`sep='\|'`, some output gene names are prefixed by taxonomy names (e.g. Pseudomonas_aeruginosa_emrE)|
|ABRicate|megares|`gene_name`|Exact match|Drop first field and `replace(':', '\|')`|
|ABRicate|[resfinder](https://bitbucket.org/genomicepidemiology/resfinder_db)|`gene_symbol`|1|`sep='_'` but gene name contains `sep`|
|[DeepARG](https://bitbucket.org/gusphdproj/deeparg-largerepo/src/master/)|deeparg|`gene_name`|Exact match||
|[ARGs-OAP](https://github.com/biofuture/Ublastx_stageone)|sarg|||[Known issue](https://github.com/AdeBC/quick_amr_db_harmonisation/issues/2)|
