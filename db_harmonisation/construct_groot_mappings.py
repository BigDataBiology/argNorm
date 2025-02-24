import pandas as pd
from argnorm.lib import get_aro_mapping_table
from Bio import SeqIO
from Bio.Seq import translate, reverse_complement, Seq
import subprocess
import requests
import os
from jug import TaskGenerator, barrier

def download_file(url, ofile):
    os.makedirs('dbs', exist_ok=True)
    os.makedirs('mapping', exist_ok=True)
    os.makedirs('manual_curation', exist_ok=True)

    with open(ofile, 'wb') as f:
        r = requests.get(url, stream=True)
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    return ofile

@TaskGenerator
def get_groot_argannot_db():
    argannot_aro_mapping_table = get_aro_mapping_table('argannot')
    argannot_aro_mapping_table.index = argannot_aro_mapping_table.index.str.lower()
            
    groot_argannot_genes = []
    groot_argannot_accessions = []
    groot_argannot_cutoffs = []
    preprocessed_argannot_refs = []

    for i in list(argannot_aro_mapping_table.index):
        gene_name = str(i).lower().split(':')
        if not str(gene_name[2][0]).isnumeric() and not '-' in gene_name[2]:
            preprocessed_argannot_refs.append(':'.join([gene_name[1], gene_name[3]]))
        else:
            preprocessed_argannot_refs.append(':'.join(i.split(':')[1:3]))
    
    argannot_aro_mapping_table.index = preprocessed_argannot_refs
    download_file('https://www.mediterranee-infection.com/wp-content/uploads/2019/03/argannot-aa-v3-march2017.txt', 'dbs/argannot_groot_db.faa')
        
    with open('./dbs/argannot_groot_db.faa', 'r') as ifile, open('./manual_curation/groot_missing.fasta', 'w') as ofile:
        for record in SeqIO.parse(ifile, 'fasta'):
            gene = str(record.id).lower().split(':')
            if not str(gene[2][0]).isnumeric() and not '-' in gene[2]:
                id = ':'.join([gene[1], gene[3]]).lower()
            else:
                id = ':'.join(str(record.id).split(':')[1:3]).lower()
                
            if id in preprocessed_argannot_refs:
                groot_argannot_genes.append(record.id)
                groot_argannot_accessions.append(str(argannot_aro_mapping_table.loc[id, 'ARO']).replace('ARO:', ''))
                groot_argannot_cutoffs.append(argannot_aro_mapping_table.loc[id, 'Cut_Off'])
            else:
                # These genes will be mapped through RGI as they aren't in the default argannot mappings
                SeqIO.write(record, ofile, 'fasta')
    argannot_groot_mapping = pd.DataFrame(groot_argannot_genes, columns=['Original ID'])
    argannot_groot_mapping['ARO'] = groot_argannot_accessions
    argannot_groot_mapping['Cut_Off'] = groot_argannot_cutoffs
    ofile = './mapping/argannot_groot_mapping.tsv'
    argannot_groot_mapping.to_csv(ofile, sep='\t', index=False)
    return ofile

@TaskGenerator
def get_groot_resfinder_db():
    os.makedirs('resfinder', exist_ok=True)
    download_file('https://bitbucket.org/genomicepidemiology/resfinder_db/get/dc33e2f9ec2c.zip', './resfinder/resfinder.zip')
    
    combine = """
        cd resfinder
        unzip resfinder.zip
        awk 'FNR==1{print ""}1' genomic*/*.fsa > resfinder_groot_db.fna
        mv resfinder_groot_db.fna ../dbs/resfinder_groot_db.fna
        cd .. && rm -rf resfinder
    """
    subprocess.run(combine, shell=True, check=True)
    
    resfinder_aro_mapping_table = get_aro_mapping_table('resfinder')
    groot_resfinder_genes = [] 
    groot_resfinder_accessions = []
    groot_resfinder_cut_offs = []

    with open('./dbs/resfinder_groot_db.fna', 'r') as ifile, open('./manual_curation/groot_missing.fasta', 'a') as ofile:
        for record in SeqIO.parse(ifile, 'fasta'):
            if str(record.id) in list(resfinder_aro_mapping_table.index):
                groot_resfinder_genes.append(record.id)
                groot_resfinder_accessions.append(str(resfinder_aro_mapping_table.loc[record.id, 'ARO']).replace('ARO:', ''))
                groot_resfinder_cut_offs.append(resfinder_aro_mapping_table.loc[record.id, 'Cut_Off'])
            else:                
                if record.seq.reverse_complement().translate().endswith('*')\
                    and record.seq.reverse_complement().translate().count('*') == 1:
                    record.seq = Seq(str(reverse_complement(translate(record.seq))).replace('*', ''))
                    SeqIO.write(record, ofile, 'fasta')
                else:
                    record.seq = Seq(str(translate(record.seq)).replace('*', ''))
                    SeqIO.write(record, ofile, 'fasta')     

    resfinder_groot_mapping = pd.DataFrame(groot_resfinder_genes, columns=['Original ID'])
    resfinder_groot_mapping['ARO'] = groot_resfinder_accessions
    resfinder_groot_mapping['Cut_Off'] = groot_resfinder_cut_offs
    
    ofile = './mapping/resfinder_groot_mapping.tsv'
    resfinder_groot_mapping.to_csv(ofile, sep='\t', index=False)
    return ofile

@TaskGenerator
def get_groot_card_db():
    os.makedirs('card', exist_ok=True)
    download_file('https://card.mcmaster.ca/latest/data', './card/data.tar.bz2')
    subprocess.check_call('cd card && tar -xvf data.tar.bz2', shell=True)
    os.rename('./card/nucleotide_fasta_protein_homolog_model.fasta', './dbs/card_groot_db.fna')
    subprocess.check_call('rm -rf card', shell=True)

    card_genes = []
    card_accessions = []
    with open('./dbs/card_groot_db.fna', 'r') as ifile:
        for record in SeqIO.parse(ifile, 'fasta'):
            card_genes.append(str(record.id).split('|')[-1])
            card_accessions.append(str(record.id).split('|')[-2].replace('ARO:', ''))
            
    card_groot_mapping = pd.DataFrame(card_genes, columns=['Original ID'])
    card_groot_mapping['ARO'] = card_accessions
    card_groot_mapping['Cut_Off'] = 'Perfect'
    
    ofile = './mapping/card_groot_mapping.tsv'
    card_groot_mapping.to_csv(ofile, sep='\t', index=False)
    return ofile

@TaskGenerator
def get_groot_missing():
    groot_missing_genes = []
    with open('./manual_curation/groot_missing.fasta', 'r') as ifile:
        for record in SeqIO.parse(ifile, 'fasta'):
            groot_missing_genes.append(record.id)

    subprocess.check_call([
        'rgi',
        'main',
        '-i', './manual_curation/groot_missing.fasta',
        '-o', './mapping/groot_missing_mapping',
        '-t', 'protein',
        '-a', 'BLAST',
        '--clean',
        '--include_loose',
        '--local'
    ])
    
    missing_groot_mappings = pd.read_csv('./mapping/groot_missing_mapping.txt', sep='\t')
    missing_groot_mappings['Original ID'] = missing_groot_mappings['ORF_ID']
    missing_groot_mappings_trim = missing_groot_mappings[['Original ID', 'ARO', 'Cut_Off']]
    ofile = './mapping/missing_groot_mapping.tsv'
    missing_groot_mappings_trim.to_csv(ofile, sep='\t', index=False)
    return ofile

@TaskGenerator
def combine_groot_mappings(argannot_path, resfinder_path, card_path, missing_path):
    argannot_groot_mapping = pd.read_csv(argannot_path, sep='\t')
    resfinder_groot_mapping = pd.read_csv(resfinder_path, sep='\t')
    card_groot_mapping = pd.read_csv(card_path, sep='\t')
    missing_groot_mapping = pd.read_csv(missing_path, sep='\t')

    comb_groot_mapping = pd.concat([
        argannot_groot_mapping,
        resfinder_groot_mapping,
        card_groot_mapping,
        missing_groot_mapping
    ]).sort_values(by=['Original ID'])
    oname_aro = 'mapping/groot_ARO_mapping.tsv'
    comb_groot_mapping.to_csv(oname_aro, sep='\t', index=False)

    groot_missing_genes = []
    with open('./manual_curation/groot_missing.fasta', 'r') as ifile:
        for record in SeqIO.parse(ifile, 'fasta'):
            groot_missing_genes.append(record.id)
            
    groot_manual_curation = pd.DataFrame(list(set(groot_missing_genes) - set(comb_groot_mapping['Original ID'])), columns=['Original ID'])
    oname_manual = 'manual_curation/groot_curation.tsv'
    groot_manual_curation.to_csv(oname_manual, sep='\t', index=False)
    return oname_aro, oname_manual

@TaskGenerator
def copy_file(oname, dest):
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    os.rename(oname, dest)
    return dest

def get_groot_aro_mapping():
    argannot_input = get_groot_argannot_db()
    resfinder_input = get_groot_resfinder_db()
    card_input = get_groot_card_db()
    missing_input = get_groot_missing()
    onames = combine_groot_mappings(argannot_input, resfinder_input, card_input, missing_input)
    copy_file(onames[0], '../argnorm/data/groot_ARO_mapping.tsv')
