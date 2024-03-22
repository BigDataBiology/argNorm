import subprocess
import pandas as pd
import os

os.makedirs('integration_tests/outputs/raw', exist_ok=True)
os.makedirs('integration_tests/outputs/hamronized', exist_ok=True)

def perform_cli_test(tool, file, folder, db=None):
    command = ['argnorm', tool, '-i', f'examples/{folder}/{file}', '-o', f'integration_tests/outputs/{folder}/{file}']
    
    if folder == 'hamronized':
        command.append('--hamronized')

    if tool == 'abricate':
        command += ['--db', db]

    subprocess.check_call(command)

    output = pd.read_csv(f'integration_tests/outputs/{folder}/{file}', sep='\t')
    golden_file = pd.read_csv(f'outputs/{folder}/{file}', sep='\t')
    pd.testing.assert_frame_equal(output, golden_file)

def test_argsoap_cli(folder):
    file = 'args-oap.sarg.reads.tsv'
    perform_cli_test('argsoap', file, folder)

test_argsoap_cli('hamronized')
test_argsoap_cli('raw')

def test_deeparg_cli(folder):
    file = 'deeparg.deeparg.orfs.tsv'
    perform_cli_test('deeparg', file, folder)

test_deeparg_cli('hamronized')
test_deeparg_cli('raw')

def test_abricate_normalizer_cli(database, folder):
    file = f'abricate.{database}.tsv'

    if not (folder == 'raw' and db == 'resfinder'):
        perform_cli_test('abricate', file, folder, db=database)

for db in ['argannot', 'megares', 'ncbi', 'resfinder']:
    test_abricate_normalizer_cli(db, 'hamronized')
    test_abricate_normalizer_cli(db, 'raw')

def test_resfinder_cli(folder, mode):
    file = f'resfinder.resfinder.{mode}.tsv'
    perform_cli_test('resfinder', file, folder)

test_resfinder_cli('hamronized', 'reads')
test_resfinder_cli('hamronized', 'orfs')
test_resfinder_cli('raw', 'reads')
test_resfinder_cli('raw', 'orfs')

def test_amrfinderplus_cli(folder):
    file = 'amrfinderplus.ncbi.orfs.tsv'
    perform_cli_test('amrfinderplus', file, folder)

test_amrfinderplus_cli('hamronized')
test_amrfinderplus_cli('raw')