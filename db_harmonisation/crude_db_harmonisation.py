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
import pandas as pd

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

def get_resfinderfg_db():
    url = 'https://raw.githubusercontent.com/RemiGSC/ResFinder_FG_Construction/606b4768433079d55f5b179219e080a45bf59dfc/output/RFG_db/ResFinder_FG.faa'
    return download_file(url, 'dbs/resfinderfg.faa')

def get_ncbi_db():
    ofile = 'dbs/ncbi_amr_raw.faa'
    url = 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.12/2024-01-31.1/AMRProt'
    return download_file(url, ofile)

def get_resfinder_db():
    ofile = 'dbs/resfinder.fna'
    url = 'https://bitbucket.org/genomicepidemiology/resfinder_db/raw/8aad1d20603fbec937cdae55024568de6dbd609f/all.fsa'
    return download_file(url, ofile)

@TaskGenerator
def get_sarg_db():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_zip = f'{tmpdir}/Short_subdatabase_V3.2.1.zip'
        download_file.f('https://smile.hku.hk/ARGs/dataset/indexingdownload/Short_subdatabase_V3.2.1.zip', tmp_zip)
        subprocess.check_call(
            ['unzip', tmp_zip, '-d', tmpdir])
        shutil.copy(f'{tmpdir}/Short_subdatabase/4.SARG_v3.2_20220917_Short_subdatabase.fasta', 'dbs/sarg.faa')

    return 'dbs/sarg.faa'

def get_deeparg_db():
    url = 'https://bitbucket.org/gusphdproj/deeparg-largerepo/raw/5683ea1c075dad3a68e0e236c98e2a98d564f560/database/v2/features.fasta'
    return download_file(url, 'dbs/deeparg.faa')

def get_argannot_db():
    url = 'https://www.mediterranee-infection.com/wp-content/uploads/2019/06/ARG_ANNOT_V5_AA_JUNE2019.txt'
    return download_file(url, 'dbs/argannot.faa')

# NCBI db has '*' at end of each protein sequence. RGI can't handle that, so '*' is removed
@TaskGenerator
def fix_ncbi(ncbi_amr_faa):
    ofile = './dbs/ncbi.faa'
    with open(ncbi_amr_faa) as original, \
            open(ofile, 'w') as corrected:
        for record in SeqIO.parse(ncbi_amr_faa, 'fasta'):
            record.seq = Seq(str(record.seq).replace("*", ""))
            SeqIO.write(record, corrected, 'fasta')

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
            '--include_loose'
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
    dfs = []

    for dir in os.listdir('./mapping'):
        if 'rgi.txt' in dir or 'rgi_output.txt' in dir:
            dfs.append(pd.read_csv(f'./mapping/{dir}', sep='\t'))
            
    comb_df = pd.concat(dfs)
    comb_df.to_csv('./mapping/combined_ARO_mapping.tsv', sep='\t')

    for i in set(comb_df['Cut_Off']):
        print(f"{i} hits: {list(comb_df['Cut_Off']).count(i) / len(list(comb_df['Cut_Off'])) * 100}% ({list(comb_df['Cut_Off']).count(i)})")

# Calling tasks
create_out_dirs()
barrier()
for db in [
        fna_to_faa(get_resfinder_db()),
        fix_ncbi(get_ncbi_db()),
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
