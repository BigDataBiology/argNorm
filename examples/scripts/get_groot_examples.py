import subprocess

# Following this documentation: https://groot-documentation.readthedocs.io/en/latest/using-groot.html#an-example

def get_input_data():
    subprocess.run("""
        echo Getting sample fastq file
        fastq-dump SRR4454613
    """, shell=True)
    
def get_groot_example(db, fastq):
    command = f"""
        echo generating GROOT {db} db example
        mkdir {db} && cd {db}
        groot get -d {db}
        groot index -m {db}.90 -i grootIndex -w 100 -p 8
        groot align -i grootIndex -f {fastq} -p 8 | groot report -c 0.95
        cd .. && rm -rf {db}
    """
    
    subprocess.check_call(command, shell=True)

get_input_data()
for db in ['resfinder', 'arg-annot', 'groot-db', 'groot-core-db', 'card']:
    get_groot_example(db, '../SRR4454613.fastq')