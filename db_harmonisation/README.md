# db_harmonisation

This is useful to _build_ the argNorm mapping databases, it is **not** intended for regular users of argNorm.

# Mapping Databases to RGI

Applicable to all databases except MEGARes

## Setting up conda environment

```bash
conda create --name rgi --channel conda-forge --channel bioconda --channel defaults rgi
conda activate rgi
conda install jug
```

## Running the code

This is a [Jug](https://jug.rtfd.io/) project, so you run `jug execute`

```bash
jug execute crude_db_harmonisation.py
```

## Constructing MEGARes Mappings

> `construct_megares_mapping.py` is integrated with `crude_db_harmonisation.py` so `crude_db_harmonisation.py` will get megares mappings as well.

MEGARes is composed of genes from:
- CARD
- NCBI
- AMRFinder
- ARG-ANNOT
- BacMet: http://bacmet.biomedicine.gu.se/, https://academic.oup.com/nar/article/42/D1/D737/1064363
- Lahey
- ResFinder
- PointFinder
- rRNACARD (same as CARD, only 2 entries have this)

```
rRNACARD,>rRNACARD|gb|NR_103957.1|+|0-2910|ARO:3004836|Ngon_23S_AZM,MEG_8254|Drugs|MLS|Macrolide-resistant_23S_rRNA_mutation|MLS23S|RequiresSNPConfirmation,

rRNACARD,>rRNACARD|gb|AE017221.1|+|1534489-1537382|ARO:3005083|Tthe_23S_PLM,MEG_8443|Drugs|Pleuromutilin|Pleuromutilin-resistant_23S_rRNA_mutation|P23S|RequiresSNPConfirmation,
```

> BacMet is ignored as metal and biocide resistance is not supported by CARD.

Megares to external header mappings file has one to many mappings!
```
- ARG-ANNOT,(Bla)SHV-1:AF148850:6-866:861,MEG_6240|Drugs|betalactams|Class_A_betalactamases|SHV
- ARG-ANNOT,(Bla)SHV-1:AF148850:6-866:861,MEG_6241|Drugs|betalactams|Class_A_betalactamases|SHV
&
- ARG-ANNOT,(Bla)OXY1-5:AY077483:181-1056:876,MEG_5284|Drugs|betalactams|Class_A_betalactamases|OXY
- ARG-ANNOT,(Bla)OXY1-5:AY077483:181-1056:876,MEG_5285|Drugs|betalactams|Class_A_betalactamases|OXY
```

Megares to external header mappings file has CARD mappings which can be directly used. However, CARD mappings have been run through RGI as there are cases where the ARO mapping for a gene from CARD is not even in the ARO:

```
CARD,CARD|gb|HM560971.1|-|119071-119269|ARO:3004468|CrpP [Pseudomonas aeruginosa] ,MEG_2133|Drugs|Fluoroquinolones|Fluoroquinolone_resistance_phosphotransferase|CRPP,

Here, ARO:3004468 does not even exist in the ARO.
```

Resfinder & argannot mappings are directly looked up using argNorm mappings. CARD, NCBI, AMRFinder, Lahey & PointFinder are ran through RGI. Mappings for resfinder & argannot that can't be found on argNorm's mappings are also run through RGI.

```
ResFinder,sitABCD_1_AY598030,MEG_8700|Multi-compound|Biocide_and_metal_resistance|Biocide_and_metal_ABC_efflux_pumps|SITABCD,

The Megares to external header mappings file states this is in resfinder, however resfinder does not have biocide & metal resistance entries and sitABCD_1_AY598030 can't be found in resfinder's mappings.
```

If a sequence translates clean or its reverse complement translates clean, the sequence is considered a CDS and the translated amino acid sequences are run through RGI using `protein` mode. Everything else is run through RGI using `contig` mode.

The RGI outputs of CDSs & contigs are combined with resfinder & argannot mappings to generate megares mappings.

> Note: as BacMet is now ignored, there will be less megares mappings (compared to when all of megares was run through RGI using `contig` mode) as loose hits from RGI won't be associated with BacMet entries.

# Groot Databases

Groot uses `resfinder`, `CARD`, and `argannot`. 
While the [groot documentation](https://groot-documentation.readthedocs.io/en/latest/groot-databases.html) says that the latest versions of card, argannot, and resfinder are used to create the groot db, only the latest version of the CARD database is downloaded, and older versions of resfinder and argannot (v3) are used - these were the latest versions of resfinder and argannot when the groot documentation was written.

- `argannot` db that groot uses: https://www.mediterranee-infection.com/wp-content/uploads/2019/03/argannot-aa-v3-march2017.txt
- `resfinder` db that groot uses: https://bitbucket.org/genomicepidemiology/resfinder_db/src/dc33e2f9ec2c420f99f77c5c33ae3faa79c999f2/

<details>
<summary>Here is the bash script that groot uses to create its database. You can view the code online <a href="https://raw.githubusercontent.com/will-rowe/groot/master/db/groot-database/make-groot-dbs.sh">here</a></summary>

```
#!/bin/env bash
#
# This script generates the groot-core-db and the groot-db
#
# It downloads the latest versions of the following databases:
#
#   - argannot
#   - resfinder
#   - card
#   - megares
#
# It then either:
#   1. groot-db:  merge the databases, remove duplicates and then cluster
#   or
#   2. groot-core-db: extract common ARGs from all databases and then cluster these common ARGs
#
# REQUIRES: Vsearch, SeqKit
#

echo "making the groot and groot-core databases..."
mkdir tmp && cd $_

# Download the latest CARD database
mkdir card && cd $_
wget --no-check-certificate -qO- https://card.mcmaster.ca/download | grep -Eo "download/0/broadstreet[a-zA-Z0-9./?=_-]*" | sort | uniq | tail -n 1 > card-db-version
cardLink=$(sed  's/^/https:\/\/card.mcmaster.ca\//g' card-db-version)
wget --no-check-certificate -O card-db.tar.gz $cardLink
tar -xvf card-db.*
awk '/>/{sub(">",">groot-db_CARD__")}1' nucleotide_fasta_protein_homolog_model.fasta > ../card-refs.fna
cd .. && rm -r card

# Download the latest ARG-annot database (V3)
wget http://en.mediterranee-infection.com/arkotheque/client/ihumed/_depot_arko/articles/1424/arg-annot-nt-v3-march2017_doc.fasta -O argannot-refs.fna
awk '/>/{sub(">",">groot-db_ARGANNOT__")}1' argannot-refs.fna > tmp && mv tmp argannot-refs.fna

# Download the latest ResFinder database
mkdir resfinder && cd $_
wget https://bitbucket.org/genomicepidemiology/resfinder_db/get/dc33e2f9ec2c.zip -O resfinder.zip
unzip resfinder.zip
awk 'FNR==1{print ""}1' genomic*/*.fsa > resfinder-refs.fna
awk '/>/{sub(">",">groot-db_RESFINDER__")}1' resfinder-refs.fna > ../resfinder-refs.fna
cd .. && rm -r resfinder

# Download the latest megres database
#mkdir megares && cd $_
#wget --no-check-certificate -qO- https://megares.meglab.org/download/index.php | grep -Eo "megares_v.*/megares_database.*[0-9].fasta" | sort | uniq | tail -n 1 > megres-db-version
#megaresLink=$(sed  's/^/https:\/\/megares.meglab.org\/download\//g' megres-db-version)
#wget --no-check-certificate -O ../megares-refs.fna $megaresLink
#cd .. && rm -r megares

# Create a reference file for the complete database
cat *.fna > all-args.fasta
seqkit rmdup --by-seq --ignore-case -j 8 -o all-args.dedup.fasta < all-args.fasta

# Cluster total set and create groot-db
mkdir groot-db.90 && cd $_
vsearch --cluster_size ../all-args.dedup.fasta --id 0.90 --msaout MSA.tmp
awk '!a[$0]++ {of="./cluster-" ++fc ".msa"; print $0 >> of ; close(of)}' RS= ORS="\n\n" MSA.tmp && rm MSA.tmp
date +%x_%H:%M:%S:%N | sed 's/\(:[0-9][0-9]\)[0-9]*$/\1/' > timestamp.txt
cd ..

# Create the groot-core-db
mkdir groot-core-db.90 && cd $_
cat ../*.fna > all-args.fasta
vsearch --cluster_size ../all-args.fasta --id 0.99 --msaout MSA.tmp
awk '!a[$0]++ {of="./cluster-" ++fc ".msa"; print $0 >> of ; close(of)}' RS= ORS="\n\n" MSA.tmp && rm MSA.tmp
for i in *.msa
do
    seqNum=$(grep '>' ${i} | wc -l)
    if [[ ${seqNum} > 3 ]]; then
        grep '>' -m 1 ${i} | sed 's/>\*//' >> core-seqs.txt
    fi
    rm ${i}
done
cat all-args.fasta | seqkit grep -f core-seqs.txt > core-seqs.fna
vsearch --cluster_size core-seqs.fna --id 0.90 --msaout MSA.tmp
awk '!a[$0]++ {of="./cluster-" ++fc ".msa"; print $0 >> of ; close(of)}' RS= ORS="\n\n" MSA.tmp && rm MSA.tmp
date +%x_%H:%M:%S:%N | sed 's/\(:[0-9][0-9]\)[0-9]*$/\1/' > timestamp.txt
rm core-seq* all-arg*
cd ..

# Finish up
mv groot* ..
cd .. && rm -r tmp
```
</details>

# Notes on the general approach

Genes from ARG annotation outputs are mapped to ARO accessions using ARO annotation tables. ARO annotation tables are constructed using the RGI (Alcock et al., 2023) and manual curation. All databases except MEGARes v3.0 are handled as amino acid files, where all ARG sequences (coding sequences) in the databases are translated to amino acid form using BioPython (Cock et al., 2009). The amino acid files are processed by RGI using the ‘protein’ mode to map genes to the ARO. The ‘Original ID’ (gene name in ARG annotation database), ‘Best_Hit_ARO’ (gene name in CARD) and ‘ARO’ (ARO accession) columns from the RGI output are specifically chosen from the RGI output, with the ‘Best_Hit_ARO’ column renamed to ‘Gene Name in CARD’ , to form the automated annotation tables. Genes which were not given an ARO mapping by RGI are manually assigned an ARO accession. The manual curation and automated annotation tables are combined to produce ARO annotation tables.

All RGI hits are checked for mapping accuracy using `check_mapping_accuracy.py`. This is done by comparing the drug categorization of the ARO assigned to the ARG with the drug categorization assigned to the ARG by its original database. A list of ARGs with mismatched mappings are included in `manual_curation/`. Mismatched hits are manually curated to correct ARGs and manual curation files are updated. Some mismatched hits are marked `correct` or `incorrect`. `correct` hits have a drug category mismatch but have an ARO mapping that's determined to be correct. `incorrect` hits have a drug category mismatch and a manually curated ARO can't be found for them. Both `correct` and `incorrect` hits should not be included in the final manual curation. For megares, a file called `megares_meta_biocide_and_virulence_genes.tsv` is also created with metal, biocide, and virulence genes. These should not be mapped to any ARO and can be manually curated to no ARO. Metal, biocide, and virulence genes are directly added to manual curation files for other databases. Groot is derived from resfinder and argannot. All loose hit manual curation for resfinder and argannot should also be manually updated for groot.

A list of all mismatched hits, with their manual curation can be found in `all_mismatched_hits.tsv`

## Handling ResFinder
The ResFinder v4.0 is a notable example as it contains forty instances of gene clusters or sequences with multiple coding sequences within (e.g. vanM_1_FJ349556 contains seven coding sequences). As RGI maps a single gene from the ResFinder database to multiple AROs (corresponding to different coding sequences within) in cases with gene clusters, gene clusters are manually assigned ARO accessions corresponding to specific gene clusters. Reverse complement sequences within ResFinder are also manually assigned ARO accessions.

## Handling MEGARes
MEGARes v3.0 is composed of multiple public genomic repositories including ResFinder, ARG-ANNOT, the National Center for Biotechnology Information (NCBI) Lahey Clinic beta-lactamase archive, the Comprehensive Antibiotic Resistance Database (CARD), NCBI’s Bacterial Antimicrobial Resistance Reference Gene Database and BacMet. The CARD database and hence the ARO do not support metal and biocide resistance, and so all entries from BacMet are removed to construct the ARO mapping table for MEGARes. Entries from ResFinder and ARG-ANNOT are directly mapped using respective ARO annotation tables, while all other entries are mapped using RGI. All coding sequences and reverse complements are translated to amino acid sequences and processed by RGI using the ‘protein’ mode. All other sequences are processed by RGI using the ‘contig’ mode. The RGI outputs are modified as mentioned above and combined with ResFinder and ARG-ANNOT mappings to form automated annotation tables. Genes which were not given an ARO mapping by RGI are manually assigned an ARO accession. The manual curation and automated annotation tables are combined to produce ARO annotation tables.

> Note: While manually assigning ARO accessions to genes, if the gene cannot be found in the ARO and if the gene does not have a parent ARO accession, no ARO mapping will be provided.

## ./mapping/megares_resfinder_argannot_mapping.tsv

Megares annotations are derived from argannot and resfinder mappings. In `../tests/test_lib.py`, `test_megares_mappings()` checks to to see if megares mappings taken directly from argannot and resfinder are correct/updated because they can be missed when running `construct_megares_mappings.py` as it depends on existing argnorm mappings.

`test_megares_mappings()` depends on `./mapping/megaresmegares_resfinder_argannot_mapping.tsv`