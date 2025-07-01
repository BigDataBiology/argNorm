from jug import Task, TaskGenerator, barrier
import shutil
import subprocess
import requests
import os
from os import path
from Bio import SeqIO
from Bio.Seq import translate, Seq
from construct_megares_mapping import construct_megares
from construct_groot_mappings import get_groot_aro_mapping


@Task
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
    subprocess.check_call(['wget', 'https://card.mcmaster.ca/download/0/broadstreet-v4.0.0.tar.bz2'])
    subprocess.check_call(['tar', '-xvf', 'broadstreet-v4.0.0.tar.bz2', './card.json'])
    subprocess.check_call(['rgi', 'load', '--card_json', 'card.json', '--local'])
    
def get_abricate_card_db():
    url = 'https://raw.githubusercontent.com/tseemann/abricate/refs/heads/master/db/card/sequences'
    return download_file(url, 'dbs/abricate_card.fna')

def get_resfinderfg_db():
    url = 'https://raw.githubusercontent.com/RemiGSC/ResFinder_FG_Construction/606b4768433079d55f5b179219e080a45bf59dfc/output/RFG_db/ResFinder_FG.faa'
    return download_file(url, 'dbs/resfinder_fg.faa')

def get_ncbi3_db():
    protein_ofile = 'dbs/ncbi3_amr_raw.faa'
    cds_ofile = 'dbs/ncbi3_amr_raw.fna'
    
    protein_url = 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.12/2024-01-31.1/AMRProt'
    cds_url = 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.12/2024-01-31.1/AMR_CDS'
    
    download_file(cds_url, cds_ofile)
    return download_file(protein_url, protein_ofile)

def get_ncbi4_db():
    protein_ofile = 'dbs/ncbi4_amr_raw.faa'
    cds_ofile = 'dbs/ncbi4_amr_raw.fna'

    protein_url = 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/4.0/2024-12-18.1/AMRProt.fa'
    cds_url = 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/4.0/2024-12-18.1/AMR_CDS.fa'
    
    download_file(cds_url, cds_ofile)
    return download_file(protein_url, protein_ofile)

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

@TaskGenerator
def combine_ncbi_dbs(ncbi3_amr_faa, ncbi4_amr_faa):
    ofile = './dbs/combined_ncbi.faa'
    with open(ncbi3_amr_faa) as ncbi3, open(ncbi4_amr_faa) as ncbi4, \
            open(ofile, 'w') as corrected:
        ncbi3_records = list(SeqIO.parse(ncbi3, 'fasta'))
        ncbi4_records = list(SeqIO.parse(ncbi4, 'fasta'))
        ncbi_records = ncbi3_records + ncbi4_records
        
        records = []
        for record in ncbi_records:
            # NCBI db has '*' at end of each protein sequence. RGI can't handle that, so '*' is removed
            record.seq = Seq(str(record.seq).replace("*", ""))
            if not record.id in records:
                SeqIO.write(record, corrected, 'fasta')
                records.append(record.id)
            else:
                continue

    return ofile

@TaskGenerator
def add_refseq_ids_to_ncbi(ncbi_faa):
    ofile = './dbs/ncbi.faa'
    ncbi3_cds = 'dbs/ncbi3_amr_raw.fna'
    ncbi4_cds = 'dbs/ncbi4_amr_raw.fna'

    with open(ncbi_faa) as ncbi_faa, open(ofile, 'w') as output:
        ncbi3_cds_records = list(SeqIO.parse(ncbi3_cds, 'fasta'))
        ncbi4_cds_records = list(SeqIO.parse(ncbi4_cds, 'fasta'))
        combined_cds_records = ncbi3_cds_records + ncbi4_cds_records
        protein_records = list(SeqIO.parse(ncbi_faa, 'fasta'))
        
        for protein_record in protein_records:
            protein_record_id = str(protein_record.id).split('|')
            for cds_record in combined_cds_records:
                if protein_record_id[1] == str(cds_record).split('|')[1]:
                    protein_record_id.insert(2, f'{str(cds_record).split('|')[2]}')
                    protein_record.id = '|'.join(protein_record_id)
                    protein_record.name = ''
                    protein_record.description = ''
                    SeqIO.write(protein_record, output, 'fasta')
                    break
            else:
                SeqIO.write(protein_record, output, 'fasta')
    
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

#
@TaskGenerator
def copy_mappings_to_argnorm(aro_mapping):
    shutil.copy(aro_mapping, '../argnorm/data')


@TaskGenerator
def get_rgi_hit_counts():
    import pandas as pd

    dfs = []
    for f in os.listdir('./mapping'):
        if 'rgi.txt' in f or 'rgi_output.txt' in f:
            dfs.append(pd.read_csv(f'./mapping/{f}', sep='\t'))

    comb_df = pd.concat(dfs)
    oname = './mapping/combined_rgi_output.tsv'
    comb_df.to_csv(oname, sep='\t', index=False)
    return oname

load_rgi()
barrier()

for db in [
        fna_to_faa(get_resfinder_db()),
        add_refseq_ids_to_ncbi(combine_ncbi_dbs(get_ncbi3_db(), get_ncbi4_db())),
        get_sarg_db(),
        get_resfinderfg_db(),
        get_deeparg_db(),
        get_argannot_db(),
        fna_to_faa(get_abricate_card_db())
    ]:
    copy_mappings_to_argnorm(run_rgi(db))
construct_megares()
get_groot_aro_mapping()
barrier()
get_rgi_hit_counts()
