from jug import TaskGenerator, barrier
import shutil
import subprocess
import requests
import os

@TaskGenerator
def create_out_dirs():
    os.makedirs('dbs', exist_ok=True)
    os.makedirs('mapping', exist_ok=True)

@TaskGenerator
def get_resfinder_db():
    from glob import glob
    subprocess.check_call(
            ['git', 'clone', 'https://bitbucket.org/genomicepidemiology/resfinder_db'])

    with open('dbs/resfinder.fna', 'w') as f:
        for file in glob('resfinder_db/*.fsa'):
            with open(file) as f2:
                f.write(f2.read())
    return 'dbs/resfinder.fna'

@TaskGenerator
def get_resfinderfg_db():
    url = 'https://raw.githubusercontent.com/RemiGSC/ResFinder_FG_Construction/606b4768433079d55f5b179219e080a45bf59dfc/output/RFG_db/ResFinder_FG.faa'
    with open('dbs/resfinder_fg.faa', 'w') as f:
        f.write(requests.get(url).text)
    return 'dbs/resfinder_fg.faa'


@TaskGenerator
def get_ncbi_db():
    subprocess.check_call(
            ['wget', 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt'])
    subprocess.check_call(
            ['mv', 'AMRProt', 'dbs/ncbi_amr.faa'])
    return 'dbs/ncbi_amr.faa'

@TaskGenerator
def get_sarg_db():
    subprocess.check_call(
            ['wget', 'https://smile.hku.hk/ARGs/dataset/indexingdownload/Short_subdatabase_V3.2.1.zip'])
    subprocess.check_call(
            ['mv', 'Short_subdatabase_V3.2.1.zip', 'sarg.zip'])
    subprocess.check_call(
            ['unzip', 'sarg.zip'])

    shutil.copy('Short_subdatabase/4.SARG_v3.2_20220917_Short_subdatabase.fasta', 'dbs/sarg.faa')
    return 'dbs/sarg.faa'

@TaskGenerator
def get_deeparg_db():
    subprocess.check_call(
            ['git', 'clone', 'https://bitbucket.org/gusphdproj/deeparg-largerepo/'])
    shutil.copy('deeparg-largerepo/database/v2/features.fasta', 'dbs/deeparg.faa')
    return 'dbs/deeparg.faa'

@TaskGenerator
def get_card_db():
    subprocess.check_call(
            ['wget', '-c', '-O', 'dbs/card.tar.bz2', 'https://card.mcmaster.ca/latest/data'])
    subprocess.check_call(
            ['tar', '-xvf', 'dbs/card.tar.bz2', '-C', 'dbs'])
    return 'dbs/card.json'

@TaskGenerator
def get_argannot_db():
    url = 'https://raw.githubusercontent.com/tseemann/abricate/master/db/argannot/sequences'
    with open('dbs/argannot.fna', 'w') as f:
        f.write(requests.get(url).text)
    return 'dbs/argannot.fna'

@TaskGenerator
def get_megares_db():
    url = 'https://www.meglab.org/downloads/megares_v3.00/megares_database_v3.00.fasta'
    with open('dbs/megares.fna', 'w') as f:
        f.write(requests.get(url).text)
    return 'dbs/megares.fna'

@TaskGenerator
def load_card_db(card_json):
    subprocess.check_call(
            ['rgi', 'load', '-i', card_json])
    
@TaskGenerator
def rgi_on_resfinder(resfinder_fna):
    subprocess.check_call(
        ['rgi', 'main', '-i', resfinder_fna, '-o', 'mapping/resfinder_rgi', '-t', 'contig', '-a', 'BLAST', '--clean', '--include_loose']
    )

@TaskGenerator
def rgi_on_ncbi(ncbi_amr_faa):
    from Bio import SeqIO
    from Bio.Seq import Seq

    with open(ncbi_amr_faa) as original, open('./dbs/ncbi_amr_corrected.faa', 'w') as corrected:
        for record in SeqIO.parse(ncbi_amr_faa, 'fasta'):
            record.seq = Seq(str(record.seq).replace("*", ""))
            SeqIO.write(record, corrected, 'fasta')

    subprocess.check_call(
        ['rgi', 'main', '-i', './dbs/ncbi_amr_corrected.faa', '-o', 'mapping/ncbi_rgi', '-t', 'protein', '-a', 'BLAST', '--clean', '--include_loose']
    )

@TaskGenerator
def rgi_on_sarg(sarg_faa):
    subprocess.check_call(
        ['rgi', 'main', '-i', sarg_faa, '-o', 'mapping/sarg_rgi', '-t', 'protein', '-a', 'BLAST', '--clean', '--include_loose']
    )

@TaskGenerator
def rgi_on_resfinderfg(resfinder_fg_faa):
    subprocess.check_call(
        ['rgi', 'main', '-i', resfinder_fg_faa, '-o', 'mapping/resfinder_fg_rgi', '-t', 'protein', '-a', 'BLAST', '--clean', '--include_loose']
    )

@TaskGenerator
def rgi_on_deeparg(deeparg_faa):
    subprocess.check_call(
        ['rgi', 'main', '-i', deeparg_faa, '-o', 'mapping/deeparg_rgi', '-t', 'protein', '-a', 'BLAST', '--clean', '--include_loose']
    )

@TaskGenerator
def rgi_on_argannot(argnnot_fna):
    subprocess.check_call(
        ['rgi', 'main', '-i', argnnot_fna, '-o', 'mapping/argannot_rgi', '-t', 'contig', '-a', 'BLAST', '--clean', '--include_loose']
    )

@TaskGenerator
def rgi_on_megares(megares_fna):
    subprocess.check_call(
        ['rgi', 'main', '-i', megares_fna, '-o', 'mapping/megares_rgi', '-t', 'contig', '-a', 'BLAST', '--clean', '--include_loose']
    )

@TaskGenerator
def reconcile_dbs():
    from .get_mapping_table import get_aro_for_hits
    
    get_aro_for_hits('mapping/resfinder_rgi.txt').to_csv('mapping/resfinder_ARO_mapping.tsv', sep='\t', index=False)
    get_aro_for_hits('mapping/ncbi_rgi.txt').to_csv('mapping/ncbi_ARO_mapping.tsv', sep='\t', index=False)
    get_aro_for_hits('mapping/sarg_rgi.txt').to_csv('mapping/sarg_ARO_mapping.tsv', sep='\t', index=False)
    get_aro_for_hits('mapping/resfinder_fg_rgi.txt').to_csv('mapping/resfinder_fg_ARO_mapping.tsv', sep='\t', index=False)
    get_aro_for_hits('mapping/deeparg_rgi.txt').to_csv('mapping/deeparg_ARO_mapping.tsv', sep='\t', index=False)
    get_aro_for_hits('mapping/argannot_rgi.txt').to_csv('mapping/argannot_ARO_mapping.tsv', sep='\t', index=False, index=False)
    get_aro_for_hits('mapping/megares_rgi.txt').to_csv('mapping/megares_ARO_mapping.tsv', sep='\t', index=False)

create_out_dirs()
databases = [
    get_resfinder_db(),
    get_resfinderfg_db(),
    get_megares_db(),
    get_ncbi_db(),
    get_sarg_db(),
    get_deeparg_db(),
    get_argannot_db()
]
card_json = get_card_db()
load_card_db(card_json)
rgi_on_databases = [
    rgi_on_resfinder('dbs/resfinder.fna'),
    rgi_on_ncbi('dbs/ncbi_amr.faa'),
    rgi_on_sarg('dbs/sarg.faa'),
    rgi_on_resfinderfg('dbs/resfinder_fg.faa'),
    rgi_on_deeparg('dbs/deeparg.faa'),
    rgi_on_argannot('dbs/argannot.fna'),
    rgi_on_megares('dbs/megares.fna')
]
reconcile_dbs()
barrier()