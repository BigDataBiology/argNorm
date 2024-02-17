# Automated AMR database harmonisation

argNorm uses two layers for AMR database harmonisation:

1. Automated annotation (implemented in this directory)
2. Manual curation


The approach used here is:

- Download NCBI and ResFinder databases
- Run them through RGI
- Using hits, create a mapping of NCBI and ResFinder genes to the ARO
- Highlight any situation of no hits/disconnect for manual curation

## Running

    conda env create -f env.yml
    conda activate crude_harmonisation
    bash crude_db_harmonisation.sh

Output mapping can be found at `resfinder_ncbi_ARO_mapping.tsv`

- 93.1% of ResFinder entries with trivial mapping to ARO using homology models (218 entries missing).
- 80.7% of NCBI entries with trivial mapping to ARO using homology models (1260 entries missing).
