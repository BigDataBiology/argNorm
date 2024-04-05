import subprocess
import pandas as pd
import os

os.makedirs('integration_tests/outputs/raw', exist_ok=True)
os.makedirs('integration_tests/outputs/hamronized', exist_ok=True)

def run_cli_test(tool, file, folder, db=None):
    command = [
        'argnorm',
        tool,
        '-i', f'examples/{folder}/{file}',
        '-o', f'integration_tests/outputs/{folder}/{file}'
    ]

    if folder == 'hamronized':
        command.append('--hamronized')

    if tool == 'abricate':
        command += ['--db', db]

    subprocess.check_call(command)

    output = pd.read_csv(f'integration_tests/outputs/{folder}/{file}', sep='\t')
    golden_file = pd.read_csv(f'outputs/{folder}/{file}', sep='\t')
    pd.testing.assert_frame_equal(output, golden_file)

for folder in ['hamronized', 'raw']:
    run_cli_test('argsoap', 'args-oap.sarg.reads.tsv', folder)
    run_cli_test('deeparg', 'deeparg.deeparg.orfs.tsv', folder)
    run_cli_test('amrfinderplus', 'amrfinderplus.ncbi.orfs.tsv', folder)
    run_cli_test('resfinder', 'resfinder.resfinder.reads.tsv', folder)
    run_cli_test('resfinder', 'resfinder.resfinder.orfs.tsv', folder)

for db in ['argannot', 'megares', 'ncbi', 'resfinder', 'resfinderfg']:
    file = f'abricate.{db}.tsv'
    run_cli_test('abricate', file, 'hamronized', db=db)
    if not 'resfinder' in db:
        run_cli_test('abricate', file, 'raw', db=db)
