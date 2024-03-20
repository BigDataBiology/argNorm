from jug import TaskGenerator, barrier
import shutil
import subprocess
import requests
import os
from os import path
import tempfile

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
    return download_file(url, 'dbs/resfinder_fg.faa')

def get_ncbi_db():
    ofile = 'dbs/ncbi_amr_raw.faa'
    url = 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt'
    return download_file(url, ofile)

@TaskGenerator
def get_resfinder_db():
    from glob import glob
    with tempfile.TemporaryDirectory() as tmpdir:
        clonedir = f'{tmpdir}/resfinder_db'
        subprocess.check_call(
            ['git', 'clone', 'https://bitbucket.org/genomicepidemiology/resfinder_db', clonedir])

        fsa_files = glob(f'{clonedir}/*.fsa')
        assert len(fsa_files) > 0, "No fsa files found"
        with open('dbs/resfinder.fna', 'w') as f:
            for file in fsa_files:
                with open(file) as f2:
                    f.write(f2.read())
    return 'dbs/resfinder.fna'


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
    url = 'https://raw.githubusercontent.com/tseemann/abricate/master/db/argannot/sequences'
    return download_file(url, 'dbs/argannot.fna')

def get_megares_db():
    url = 'https://www.meglab.org/downloads/megares_v3.00/megares_database_v3.00.fasta'
    return download_file(url, 'dbs/megares.fna')

@TaskGenerator
def fix_ncbi(ncbi_amr_faa):
    from Bio import SeqIO
    from Bio.Seq import Seq

    ofile = './dbs/ncbi.faa'
    with open(ncbi_amr_faa) as original, \
            open(ofile, 'w') as corrected:
        for record in SeqIO.parse(ncbi_amr_faa, 'fasta'):
            record.seq = Seq(str(record.seq).replace("*", ""))
            SeqIO.write(record, corrected, 'fasta')
    return ofile

@TaskGenerator
def run_rgi(fa):
    from get_mapping_table import get_aro_for_hits
    db = path.basename(fa).split('.')[0]
    if fa.endswith('.fna'):
        mode = 'contig'
    elif fa.endswith('.faa'):
        mode = 'protein'
    else:
        raise ValueError(f"Unknown file type {fa}")
    rgi_ofile = f'mapping/{db}_rgi'
    ofile = f'mapping/{db}_ARO_mapping.tsv'
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
    get_aro_for_hits(rgi_ofile + '.txt', db).to_csv(ofile, sep='\t', index=False)
    return ofile

@TaskGenerator
def move_mappings_to_argnorm(aro_mapping):
    shutil.copy(aro_mapping, '../argnorm/data')

create_out_dirs()
barrier()
for db in [
        get_resfinder_db(),
        fix_ncbi(get_ncbi_db()),
        get_sarg_db(),
        get_resfinderfg_db(),
        get_deeparg_db(),

        get_megares_db(),
        get_argannot_db()
    ]:
    move_mappings_to_argnorm(run_rgi(db))