import pandas as pd
from argnorm.lib import get_aro_mapping_table
from Bio import SeqIO
from Bio.Seq import translate
import subprocess
import requests
import os

def download_file(url, ofile):
    os.makedirs('dbs', exist_ok=True)
    os.makedirs('mapping', exist_ok=True)

    with open(ofile, 'wb') as f:
        r = requests.get(url, stream=True)
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    return ofile

argannot_aro_mapping_table = get_aro_mapping_table('argannot')
        
groot_argannot_genes = []
groot_argannot_accessions = []

argannot_aro_mapping_table.index = argannot_aro_mapping_table.index.str.lower()
replacement_aro_index = []

default_argannot_mappings = []
for i in list(argannot_aro_mapping_table.index):
    gene = str(i).lower().split(':')
    if not str(gene[2][0]).isnumeric() and not '-' in gene[2]:
        default_argannot_mappings.append(':'.join([gene[1], gene[3]]))
    else:
        default_argannot_mappings.append(':'.join(i.split(':')[1:3]))
        
argannot_aro_mapping_table.index = default_argannot_mappings

download_file('https://www.mediterranee-infection.com/wp-content/uploads/2019/03/argannot-aa-v3-march2017.txt', 'dbs/argannot_groot_db.faa')
    
with open('./dbs/argannot_groot_db.faa', 'r') as ifile, open('./groot_missing.fasta', 'w') as ofile:
    for record in SeqIO.parse(ifile, 'fasta'):
        gene = str(record.id).lower().split(':')
        if not str(gene[2][0]).isnumeric() and not '-' in gene[2]:
            id = ':'.join([gene[1], gene[3]]).lower()
        else:
            id = ':'.join(str(record.id).split(':')[1:3]).lower()
            
        if id in default_argannot_mappings:
            groot_argannot_genes.append(record.id)
            groot_argannot_accessions.append(str(argannot_aro_mapping_table.loc[id, 'ARO']).removeprefix('ARO:'))
        else:
            # Run this through RGI
            SeqIO.write(record, ofile, 'fasta')

argannot_groot_mapping = pd.DataFrame(groot_argannot_genes, columns=['Original ID'])
argannot_groot_mapping['ARO'] = groot_argannot_accessions

subprocess.check_call(['bash', 'get_groot_dbs.sh'])

resfinder_aro_mapping_table = get_aro_mapping_table('resfinder')
groot_resfinder_genes = [] 
groot_resfinder_accessions = []

with open('./resfinder-refs.fna', 'r') as ifile, open('./groot_missing.fasta', 'a') as ofile:
    for record in SeqIO.parse(ifile, 'fasta'):
        if str(record.id) in list(resfinder_aro_mapping_table.index):
            groot_resfinder_genes.append(record.id)
            groot_resfinder_accessions.append(str(resfinder_aro_mapping_table.loc[record.id, 'ARO']).removeprefix('ARO:'))
        else:
            record.seq = translate(record.seq).removesuffix('*')
            print(record)
            SeqIO.write(record, ofile, 'fasta')
            
resfinder_groot_mapping = pd.DataFrame(groot_resfinder_genes, columns=['Original ID'])
resfinder_groot_mapping['ARO'] = groot_resfinder_accessions

card_genes = [] 
card_accessions = []
with open('./card-refs.fna', 'r') as ifile:
    for record in SeqIO.parse(ifile, 'fasta'):
        card_genes.append(str(record.id).split('|')[-1])
        card_accessions.append(str(record.id).split('|')[-2].removeprefix('ARO:'))
        
card_groot_mapping = pd.DataFrame(card_genes, columns=['Original ID'])
card_groot_mapping['ARO'] = card_accessions

# subprocess.check_call([
#     'rgi',
#     'main',
#     '-i', './groot_missing.fasta',
#     '-o', './mapping/groot_missing_mapping',
#     '-t', 'protein',
#     '-a', 'BLAST',
#     '--clean',
#     '--include_loose'
# ])

missing_groot_mappings = pd.read_csv('./mapping/groot_missing_mapping.txt', sep='\t')
missing_groot_mappings['Original ID'] = missing_groot_mappings['ORF_ID']
missing_groot_mappings_trim = missing_groot_mappings[['Original ID', 'ARO']]

comb_groot_mapping = pd.concat([argannot_groot_mapping, resfinder_groot_mapping, card_groot_mapping, missing_groot_mappings_trim])
comb_groot_mapping.to_csv('./mapping/groot_ARO_mapping.tsv', sep='\t', index=False)