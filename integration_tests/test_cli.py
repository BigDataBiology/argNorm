import subprocess
import pandas as pd
import os

os.makedirs('integration_tests/outputs/raw', exist_ok=True)
os.makedirs('integration_tests/outputs/hamronized', exist_ok=True)

def test_cli(tool, file, folder, db=None):
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
    test_cli('argsoap', 'args-oap.sarg.reads.tsv', folder)
    test_cli('deeparg', 'deeparg.deeparg.orfs.tsv', folder)
    test_cli('amrfinderplus', 'amrfinderplus.ncbi.orfs.tsv', folder)
    test_cli('resfinder', 'resfinder.resfinder.reads.tsv', folder)
    test_cli('resfinder', 'resfinder.resfinder.orfs.tsv', folder)

for db in ['argannot', 'megares', 'ncbi', 'resfinder']:
    file = f'abricate.{db}.tsv'
    test_cli('abricate', file, 'hamronized', db=db)
    if db != 'resfinder':
        test_cli('abricate', file, 'raw', db=db)