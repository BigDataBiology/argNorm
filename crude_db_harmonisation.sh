#!/bin/bash

# set up dirs
mkdir -p dbs

# get the resfinder database (skipping pointfinder for now)
git clone https://bitbucket.org/genomicepidemiology/resfinder_db
cat resfinder_db/*.fsa > dbs/resfinder.fna
rm -rf resfinder_db

# get the NCBI databases
wget https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt
mv AMRProt dbs/ncbi_amr.faa

# get latest CARD database
wget -O dbs/card.tar.bz2 https://card.mcmaster.ca/latest/data 
tar -xvf dbs/card.tar.bz2 -C dbs

# load latest card into rgi
rgi load -i dbs/card.json

# run rgi on both of these databases
# got to do CDS as resfinder doesn't have 
mkdir -p mapping
rgi main -i dbs/resfinder.fna -o mapping/resfinder_rgi -t contig -a BLAST --clean
rgi main -i dbs/ncbi_amr.faa -o mapping/ncbi_rgi -t protein -a BLAST --clean

# reconcile the databases
python reconcile.py -f dbs/resfinder.fna -r mapping/resfinder_rgi.txt -d resfinder
python reconcile.py -f dbs/ncbi_amr.faa -r mapping/ncbi_rgi.txt -d ncbi

# combine outputs
awk 'FNR > 1' resfinder_ARO_mapping.tsv ncbi_ARO_mapping.tsv > resfinder_ncbi_ARO_mapping.tsv

# tidy up
mv resfinder_ARO_mapping.tsv ncbi_ARO_mapping.tsv mapping
