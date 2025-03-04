from jug import TaskGenerator, barrier
import shutil
import subprocess
import requests
import os
from os import path
import tempfile
from Bio import SeqIO
from Bio.Seq import translate, Seq
from construct_megares_mapping import construct_megares
from construct_groot_mappings import get_groot_aro_mapping

@TaskGenerator
def create_out_dirs():
    os.makedirs('dbs', exist_ok=True)
    os.makedirs('mapping', exist_ok=True)

@TaskGenerator
def download_file(url, ofile):
    with open(ofile, 'wb') as f:
        r = requests.get(url, stream=True)
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    return ofile

@TaskGenerator
def load_rgi():
    subprocess.check_call(['rgi', 'clean', '--local'])
    subprocess.check_call(['wget', 'https://card.mcmaster.ca/latest/data'])
    subprocess.check_call(['tar', '-xvf', 'data', './card.json'])
    subprocess.check_call(['rgi', 'load', '--card_json', 'card.json', '--local'])

def get_resfinderfg_db():
    url = 'https://raw.githubusercontent.com/RemiGSC/ResFinder_FG_Construction/606b4768433079d55f5b179219e080a45bf59dfc/output/RFG_db/ResFinder_FG.faa'
    return download_file(url, 'dbs/resfinder_fg.faa')

def get_ncbi3_db():
    ofile = 'dbs/ncbi3_amr_raw.faa'
    url = 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.12/2024-01-31.1/AMRProt'
    return download_file(url, ofile)

def get_ncbi4_db():
    ofile = 'dbs/ncbi4_amr_raw.faa'
    url = 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/4.0/2024-12-18.1/AMRProt.fa'
    return download_file(url, ofile)

def get_resfinder_db():
    ofile = 'dbs/resfinder.fna'
    url = 'https://bitbucket.org/genomicepidemiology/resfinder_db/raw/8aad1d20603fbec937cdae55024568de6dbd609f/all.fsa'
    return download_file(url, ofile)

def get_sarg_db():
    url = 'https://raw.githubusercontent.com/xinehc/args_oap/a3e5cff4a6c09f81e4834cfd9a31e6ce7d678d71/src/args_oap/db/sarg.fasta'
    return download_file(url, 'dbs/sarg.faa')

def get_deeparg_db():
    url = 'https://bitbucket.org/gusphdproj/deeparg-largerepo/raw/5683ea1c075dad3a68e0e236c98e2a98d564f560/database/v2/features.fasta'
    return download_file(url, 'dbs/deeparg.faa')

def get_argannot_db():
    url = 'https://www.mediterranee-infection.com/wp-content/uploads/2019/06/ARG_ANNOT_V5_AA_JUNE2019.txt'
    return download_file(url, 'dbs/argannot.faa')

# NCBI db has '*' at end of each protein sequence. RGI can't handle that, so '*' is removed
@TaskGenerator
def fix_ncbi(ncbi3_amr_faa, ncbi4_amr_faa):
    ofile = './dbs/ncbi.faa'
    with open(ncbi3_amr_faa) as ncbi3, open(ncbi4_amr_faa) as ncbi4, \
            open(ofile, 'w') as corrected:
        ncbi3_records = list(SeqIO.parse(ncbi3, 'fasta'))
        ncbi4_records = list(SeqIO.parse(ncbi4, 'fasta'))
        ncbi_records = ncbi3_records + ncbi4_records
        
        records = []
        for record in ncbi_records:
            record.seq = Seq(str(record.seq).replace("*", ""))
            if not record.id in records:
                SeqIO.write(record, corrected, 'fasta')
                records.append(record.id)
            else:
                continue

    return ofile

# Needed when nucleotide database (eg. resfinder) needs to be run through RGI with protein mode
@TaskGenerator
def fna_to_faa(ifile):
    ofile = ifile.replace('.fna', '.faa')
    with open(ifile) as original, open(ofile, 'w') as output:
        for record in SeqIO.parse(original, 'fasta'):
            record.seq = Seq(str(translate(record.seq)).replace('*', ''))
            SeqIO.write(record, output, 'fasta')

    return ofile

@TaskGenerator
def run_rgi(fa):
    from get_mapping_table import get_aro_for_hits

    db_name = path.basename(fa).split('.')[0]

    if fa.endswith('.fna'):
        mode = 'contig'
    elif fa.endswith('.faa'):
        mode = 'protein'
    else:
        raise ValueError(f"Unknown file type {fa}")

    rgi_ofile = f'mapping/{db_name}_rgi'
    ofile = f'mapping/{db_name}_ARO_mapping.tsv'

    subprocess.check_call(
        [
            'rgi',
            'main',
            '-i', fa,
            '-o', rgi_ofile,
            '-t', mode,
            '-a', 'BLAST',
            '--clean',
            '--include_loose',
            '--local'
        ]
    )

    get_aro_for_hits(fa, rgi_ofile + '.txt', db_name).to_csv(ofile, sep='\t', index=False)
    return ofile

# Moving ARO mapping tables over to argnorm/data
@TaskGenerator
def move_mappings_to_argnorm(aro_mapping):
    shutil.copy(aro_mapping, '../argnorm/data')

@TaskGenerator
def get_rgi_hit_counts():
    import pandas as pd

    dfs = []
    for dir in os.listdir('./mapping'):
        if 'rgi.txt' in dir or 'rgi_output.txt' in dir:
            dfs.append(pd.read_csv(f'./mapping/{dir}', sep='\t'))
            
    comb_df = pd.concat(dfs)
    comb_df.to_csv('./mapping/combined_ARO_mapping.tsv', sep='\t')

# Calling tasks
create_out_dirs()
load_rgi()
barrier()
for db in [
        fna_to_faa(get_resfinder_db()),
        fix_ncbi(get_ncbi3_db(), get_ncbi4_db()),
        get_sarg_db(),
        get_resfinderfg_db(),
        get_deeparg_db(),
        get_argannot_db()
    ]:
    move_mappings_to_argnorm(run_rgi(db))
construct_megares()
get_groot_aro_mapping()
barrier()
get_rgi_hit_counts()
