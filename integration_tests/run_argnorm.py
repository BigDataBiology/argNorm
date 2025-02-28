import subprocess
import pandas as pd
import os
from resource import *

os.makedirs('integration_tests/outputs/raw', exist_ok=True)
os.makedirs('integration_tests/outputs/hamronized', exist_ok=True)

def run_cli_test(tool, file, folder, db=None):
    if folder == 'hamronized':
        tool = 'hamronization'

    command = [
        'argnorm',
        tool,
        '-i', f'examples/{folder}/{file}',
        '-o', f'integration_tests/outputs/{folder}/{file}'
    ]

    if tool in  ['abricate', 'groot']:
        command += ['--db', db]
    
    baseline = getrusage(RUSAGE_SELF).ru_maxrss
    subprocess.check_call(command)
    print(f'Time taken to run for {tool} {file} {folder} : {getrusage(RUSAGE_SELF).ru_utime}')
    print(f'Baseline memory: {baseline}')
    print(f'Program memory for {tool} {file} {folder} : {getrusage(RUSAGE_SELF).ru_maxrss - baseline}')

    output = pd.read_csv(f'integration_tests/outputs/{folder}/{file}', sep='\t')
    golden_file = pd.read_csv(f'outputs/{folder}/{file}', sep='\t')
    pd.testing.assert_frame_equal(output, golden_file)

for folder in ['hamronized', 'raw']:
    run_cli_test('ARGSOAP', 'args-oap.sarg.reads.tsv', folder)
    run_cli_test('argsoap', 'args-oap.sarg.reads.tsv', folder)
    run_cli_test('DEEParg', 'deeparg.deeparg.orfs.tsv', folder)
    run_cli_test('deeparg', 'deeparg.deeparg.orfs.tsv', folder)
    run_cli_test('amrfinderplus', 'amrfinderplus.ncbi.orfs.tsv', folder)
    run_cli_test('resfinder', 'resfinder.resfinder.reads.tsv', folder)
    run_cli_test('resfinder', 'resfinder.resfinder.orfs.tsv', folder)
    for db in ['argannot', 'card', 'groot-core-db', 'groot-db', 'resfinder']:
        database = db
        if not 'groot' in db:
            database = f'groot-{db}'
        run_cli_test('groot', f'groot.{db}.tsv', folder, database)

for db in ['ARGANNOT', 'argannot', 'MEGAres', 'megares', 'ncbi', 'resfinder', 'resfinderfg']:
    file  = f'abricate.{db.lower()}.tsv'
    run_cli_test('hamronization', file, 'hamronized', db=db)
    if not 'resfinder' in db:
        run_cli_test('abricate', file, 'raw', db=db)

run_cli_test('hamronization', 'combined_hamronization.tsv', 'hamronized')