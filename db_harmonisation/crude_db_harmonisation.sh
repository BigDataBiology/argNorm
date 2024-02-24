#!/bin/bash

# Set up database directory
mkdir -p dbs

# Get resfinder db
git clone https://bitbucket.org/genomicepidemiology/resfinder_db
cat resfinder_db/*.fsa > dbs/resfinder.fna
rm -rf resfinder_db

# Get resfinderfg db
git clone https://github.com/BigDataBiology/ResFinderFG
cp ResFinderFG/ResFinder_FG.txt dbs/resfinder_fg.faa

# Get NCBI db
wget https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt
mv AMRProt dbs/ncbi_amr.faa

# Get the SARG db
wget https://smile.hku.hk/ARGs/dataset/indexingdownload/Short_subdatabase_V3.2.1.zip
mv Short_subdatabase_V3.2.1.zip sarg.zip
unzip sarg.zip
cat Short_subdatabase/4.SARG_v3.2_20220917_Short_subdatabase.fasta > dbs/sarg.faa
rm -rf Short_subdatabase
rm -rf __MACOSX
rm -rf sarg.zip

# Get DeepARG db
git clone https://bitbucket.org/gusphdproj/deeparg-largerepo/
cp deeparg-largerepo/database/v2/features.fasta dbs/deeparg.faa

# Get CARD db
wget -c -O dbs/card.tar.bz2 https://card.mcmaster.ca/latest/data
tar -xvf dbs/card.tar.bz2 -C dbs
ls -hl dbs

# Get megares db
wget -c -O dbs/megares.fna https://megares.meglab.org/download/megares_v2.00/megares_full_database_v2.00.fasta

# Get argannot db
wget -c -O dbs/argannot.fna https://raw.githubusercontent.com/tseemann/abricate/master/db/argannot/sequences

# Load latest card into rgi
rgi load -i dbs/card.json

# Running rgi on dbs
mkdir -p mapping

rgi main -i dbs/resfinder.fna -o mapping/resfinder_rgi -t contig -a BLAST --clean --include_loose
rgi main -i dbs/ncbi_amr.faa -o mapping/ncbi_rgi -t protein -a BLAST --clean --include_loose
rgi main -i dbs/sarg.faa -o mapping/sarg_rgi -t protein -a BLAST --clean --include_loose
rgi main -i dbs/resfinder_fg.faa -o mapping/resfinder_fg_rgi -t protein -a BLAST --clean --include_loose
rgi main -i dbs/deeparg.faa -o mapping/deeparg_rgi -t protein -a BLAST --clean --include_loose
rgi main -i dbs/argannot.fna -o mapping/argannot_rgi -t contig -a BLAST --clean --include_loose
rgi main -i dbs/megares.fna -o mapping/megares_rgi -t contig -a BLAST --clean --include_loose

# Reconcile the databases
python reconcile.py -f dbs/resfinder.fna -r mapping/resfinder_rgi.txt -d resfinder
python reconcile.py -f dbs/ncbi_amr.faa -r mapping/ncbi_rgi.txt -d ncbi
python reconcile.py -f dbs/sarg.faa -r mapping/sarg_rgi.txt -d sarg
python reconcile.py -f dbs/resfinder_fg.faa -r mapping/resfinder_fg_rgi.txt -d resfinder_fg
python reconcile.py -f dbs/deeparg.faa -r mapping/deeparg_rgi.txt -d deeparg
python reconcile.py -f dbs/argannot.fna -r mapping/argannot_rgi.txt -d argannot
python reconcile.py -f dbs/megares.fna -r mapping/megares_rgi.txt -d megares
