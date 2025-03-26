import subprocess
import pandas as pd
import os

os.makedirs('integration_tests/outputs/raw', exist_ok=True)
os.makedirs('integration_tests/outputs/hamronized', exist_ok=True)

def run_cli_test(tool, file, folder, *, db=None, options=[]):
    if folder == 'hamronized':
        tool = 'hamronization'

    command = [
        '/usr/bin/time',
            '-f', "Time %e seconds (%M KB mem used)",  # %e is wall clock time, %M is max memory
        'argnorm',
        tool,
        '-i', f'examples/{folder}/{file}',
        '-o', f'integration_tests/outputs/{folder}/{file}'
    ]

    if tool in ['abricate', 'groot']:
        command += ['--db', db]

    command += options
    subprocess.check_call(command)

    output = pd.read_csv(f'integration_tests/outputs/{folder}/{file}', sep='\t')
    golden_file = pd.read_csv(f'outputs/{folder}/{file}', sep='\t')
    pd.testing.assert_frame_equal(output, golden_file)

print('Testing argnorm --version')
subprocess.check_call(['argnorm', '--version'])

for folder in ['hamronized', 'raw']:
    run_cli_test('ARGSOAP', 'args-oap.sarg.reads.tsv', folder)
    run_cli_test('argsoap', 'args-oap.sarg.reads.tsv', folder)
    run_cli_test('DEEParg', 'deeparg.deeparg.orfs.tsv', folder)
    run_cli_test('deeparg', 'deeparg.deeparg.orfs.tsv', folder)
    run_cli_test('amrfinderplus', 'amrfinderplus.ncbi.orfs.v3.10.30.tsv', folder)
    run_cli_test('amrfinderplus', 'amrfinderplus.ncbi.orfs.v4.tsv', folder)
    run_cli_test('resfinder', 'resfinder.resfinder.reads.tsv', folder)
    run_cli_test('resfinder', 'resfinder.resfinder.orfs.tsv', folder)
    for db in ['argannot', 'card', 'groot-core-db', 'groot-db', 'resfinder']:
        database = db
        if not 'groot' in db:
            database = f'groot-{db}'
        run_cli_test('groot', f'groot.{db}.tsv', folder, db=database)

for db in ['ARGANNOT', 'argannot', 'MEGAres', 'megares', 'ncbi', 'resfinder', 'resfinderfg']:
    file  = f'abricate.{db.lower()}.tsv'
    run_cli_test('hamronization', file, 'hamronized', db=db)
    if not 'resfinder' in db:
        run_cli_test('abricate', file, 'raw', db=db)

run_cli_test('hamronization', 'combined_hamronization.tsv', 'hamronized')
run_cli_test('hamronization', 'combined_hamronization_full.tsv', 'hamronized', options=['--hamronization_skip_unsupported_tool'])
run_cli_test('hamronization', 'abricate.ncbi.hamronizationv4.tsv', 'hamronized')

try:
    run_cli_test('hamronization', 'combined_hamronization_full.tsv', 'hamronized')
except:
    pass
else:
    raise AssertionError(f'combined_hamronization_full.tsv is normalized without exception without using --hamronization_skip_unsupported_tool')
