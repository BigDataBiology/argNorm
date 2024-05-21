import pandas as pd
from argnorm.lib import get_aro_mapping_table
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import tempfile
import shutil
import requests
import zipfile
import os
from jug import TaskGenerator, barrier

def download_file(url, ofile):
    os.makedirs('dbs', exist_ok=True)
    os.makedirs('mapping', exist_ok=True)

    with open(ofile, 'wb') as f:
        r = requests.get(url, stream=True)
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    return ofile

def get_megares_db():
    """
    Gets the megares database and megares_to_external_header_mappings csv file which
    gives data on which databases each gene was taken from.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_zip = f'{tmpdir}/megares.zip'
        download_file('https://www.meglab.org/downloads/megares_v3.00.zip', tmp_zip)

        with zipfile.ZipFile(tmp_zip, 'r') as zip_ref:
            zip_ref.extractall(f'{tmpdir}/megares/')

        shutil.copy(f'{tmpdir}/megares/megares_v3.00/megares_database_v3.00.fasta', 'dbs/megares.fasta')
        shutil.copy(f'{tmpdir}/megares/megares_v3.00/megares_to_external_header_mappings_v3.00.csv', 'dbs/megares_headers.csv')

    return ['dbs/megares.fasta', 'dbs/megares_headers.csv']

def search_argnorm_mappings(mappings, db):
    mappings.drop(columns=['UpdatedHeader', 'Database'], inplace=True)
    aros = []
    mapping_table = get_aro_mapping_table(db)
    mapping_table.index = list(mapping_table.index.map(lambda x: str(x).lower()))

    for header in mappings['Source_header']:
        try:
            aros.append(str(mapping_table.loc[header.strip().lower(), 'ARO'])[4:])
        except KeyError:
            aros.append(np.nan)

    mappings['ARO'] = aros
    mappings.rename(columns={'MEGARes_header': 'Original ID'}, inplace=True)
    return mappings

def get_resfinder_mappings(megares_headers):
    mappings = search_argnorm_mappings(megares_headers[
        ((megares_headers.Database == 'ResFinder') | (megares_headers.Database == 'Resfinder'))\
        & (~ megares_headers.MEGARes_header.str.contains('Multi-compound'))\
        & (~ megares_headers.MEGARes_header.str.contains('Metals'))\
        & (~ megares_headers.MEGARes_header.str.contains('Biocides'))], 'resfinder')
    return mappings

def get_argannot_mappings(megares_headers):
    mappings = search_argnorm_mappings(megares_headers[megares_headers.Database == 'ARG-ANNOT'], 'argannot')
    return mappings

def get_missing_mappings(resfinder_argannot_mappings, megares_headers):
    missing_mappings = []
    mapping_table = pd.concat([get_aro_mapping_table('resfinder'), get_aro_mapping_table('argannot')])
    mapping_table.index = list(mapping_table.index.map(lambda x: str(x).lower()))

    for header in resfinder_argannot_mappings['Source_header']:
        try:
            aro = mapping_table.loc[header.strip().lower(), 'ARO']
        except KeyError:
            if type(resfinder_argannot_mappings.set_index('Source_header').loc[header, 'Original ID']) == str:
                missing_mappings.append(resfinder_argannot_mappings.set_index('Source_header').loc[header, 'Original ID'])
            else:
                missing_mappings.append(resfinder_argannot_mappings.set_index('Source_header').loc[header, 'Original ID'][0])

    # Update missing mappings to include entries from NCBI, AMRFinder, PointFinder, Lahey
    missing_mappings += list(megares_headers[(megares_headers.Database == 'PointFinder')\
                                        | (megares_headers.Database == 'AMRFinder')\
                                        | (megares_headers.Database == 'NCBI')\
                                        | (megares_headers.Database == 'Lahey')\
                                        | (megares_headers.Database == 'CARD')\
                                        | (megares_headers.Database == 'rRNACARD')]['MEGARes_header'])
    return list(set(missing_mappings))

def generate_missing_mappings_fasta(missing_mappings, megares_fasta):
    """
    Get protein file for missing CDSs to pass to RGI and nucleotide file for missing contigs to pass to RGI
    """
    ofile1 = './megares_cds.fasta'
    ofile2 = './megares_contigs.fasta'
    with open(megares_fasta) as ifile, open(ofile1, 'w') as megares_cds, open(ofile2, 'w') as megares_contigs:
        for record in SeqIO.parse(ifile, 'fasta'):
            if record.id not in missing_mappings:
                continue

            if record.seq.translate().endswith('*') and record.seq.translate().count('*') == 1:
                record.seq = Seq(str(record.seq.translate()).replace('*', ''))
                SeqIO.write(record, megares_cds, 'fasta')
            elif record.seq.reverse_complement().translate().endswith('*')\
                and record.seq.reverse_complement().translate().count('*') == 1:
                record.seq = Seq(str(record.seq.reverse_complement().translate()).replace('*', ''))
                SeqIO.write(record, megares_cds, 'fasta')
            else:
                SeqIO.write(record, megares_contigs, 'fasta')
    return [ofile1, ofile2]

@TaskGenerator
def setup_for_rgi():
    get_megares_db()
    megares_headers = pd.read_csv('./dbs/megares_headers.csv')

    resfinder_argannot_mapping = pd.concat([
        get_resfinder_mappings(megares_headers),
        get_argannot_mappings(megares_headers)
    ])
    missing_mappings = get_missing_mappings(resfinder_argannot_mapping, megares_headers)
    resfinder_argannot_mapping.drop(columns=['Source_header'], inplace=True)
    resfinder_argannot_mapping.to_csv('./mapping/megares_card_resfinder_argannot_mapping.tsv', sep='\t', index=False)
    return generate_missing_mappings_fasta(missing_mappings, './dbs/megares.fasta')

@TaskGenerator
def get_cds_rgi_output(cds_fasta):
    ofile = './megares_cds_rgi_output'
    subprocess.check_call([
        'rgi',
        'main',
        '-i', cds_fasta,
        '-o', ofile,
        '-t', 'protein',
        '-a', 'BLAST',
        '--clean',
        '--include_loose'
    ])
    return f'{ofile}.txt'

@TaskGenerator
def get_contig_rgi_output(contig_fasta):
    ofile = './megares_contigs_rgi_output'
    subprocess.check_call([
        'rgi',
        'main',
        '-i', contig_fasta,
        '-o', ofile,
        '-t', 'contig',
        '-a', 'BLAST',
        '--clean',
        '--include_loose'
    ])
    return f'{ofile}.txt'

@TaskGenerator
def merge_megares_mappings(cds_mapping, contig_mapping):
    megares_mappings = pd.read_csv('./mapping/megares_card_resfinder_argannot_mapping.tsv', sep='\t')
    cds_mapping = pd.read_csv(cds_mapping, sep='\t')
    contig_mapping = pd.read_csv(contig_mapping, sep='\t')

    print(cds_mapping)
    print(contig_mapping)

    cds_mapping['Original ID'] = cds_mapping['ORF_ID']
    contig_mapping['Original ID'] = contig_mapping['Contig'].apply(lambda x: '_'.join(x.split('_')[:-1]))

    cds_mapping = cds_mapping[['Original ID', 'ARO']]
    contig_mapping = contig_mapping[['Original ID', 'ARO']]

    print(megares_mappings)

    megares_mappings.set_index('Original ID').drop(index=set(cds_mapping['Original ID']) & set(megares_mappings.index), inplace=True)
    megares_mappings.set_index('Original ID').drop(index=set(contig_mapping['Original ID']) & set(megares_mappings.index), inplace=True)
    megares_mappings = pd.concat([megares_mappings, cds_mapping, contig_mapping])
    megares_mappings = megares_mappings.dropna().drop_duplicates()
    megares_mappings = megares_mappings.astype({'ARO': 'str'})
    megares_mappings['ARO'] = megares_mappings['ARO'].apply(lambda x: str(x).replace('.0', ''))
    megares_mappings['Database'] = 'megares'
    megares_mappings.to_csv('./mapping/megares_ARO_mapping.tsv', sep='\t', index=False)
    shutil.copy('./mapping/megares_ARO_mapping.tsv', '../argnorm/data/megares_ARO_mapping.tsv')
    return megares_mappings

@TaskGenerator
def get_megares_manual_curation(megares_mappings):
    os.makedirs('manual_curation', exist_ok=True)
    megares_headers = pd.read_csv('./dbs/megares_headers.csv')
    manual_curation = []
    megares_headers = megares_headers[~(megares_headers.Database == 'BacMet')]
    manual_curation += list(set(megares_headers['MEGARes_header']) - set(megares_mappings['Original ID']))

    one_to_many = []
    for i in list(megares_mappings['Original ID']):
        if list(megares_mappings['Original ID']).count(i) > 1:
            one_to_many.append(i)

    manual_curation += list(set(one_to_many))
    pd.DataFrame({'Original ID': manual_curation}).to_csv('./manual_curation/megares_curation.tsv', sep='\t', index=False)

def construct_megares():
    fasta_files = setup_for_rgi()
    cds_mapping = get_cds_rgi_output(fasta_files[0])
    contig_mapping = get_contig_rgi_output(fasta_files[1])
    megares_mapping = merge_megares_mappings(cds_mapping, contig_mapping)
    get_megares_manual_curation(megares_mapping)