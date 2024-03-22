import subprocess
import pytest
import pandas as pd

def perform_cli_test(tool, file, folder, db=None):
    command = ['argnorm', tool, '-i', f'examples/{folder}/{file}', '-o', f'tests/outputs/{folder}/{file}']

    if folder == 'hamronized':
        command.append('--hamronized')

    if tool == 'abricate':
        command += ['--db', db]

    subprocess.check_call(command)

    output = pd.read_csv(f'tests/outputs/{folder}/{file}', sep='\t')
    golden_file = pd.read_csv(f'outputs/{folder}/{file}', sep='\t')
    pd.testing.assert_frame_equal(output, golden_file)

@pytest.mark.parametrize('folder', ['hamronized', 'raw'])
def test_argsoap_cli(folder):
    file = 'args-oap.sarg.reads.tsv'
    perform_cli_test('argsoap', file, folder)

@pytest.mark.parametrize('folder', ['hamronized', 'raw'])
def test_deeparg_cli(folder):
    file = 'deeparg.deeparg.orfs.tsv'
    perform_cli_test('deeparg', file, folder)

@pytest.mark.parametrize('database', ['argannot', 'megares', 'ncbi', 'resfinder'])
def test_abricate_normalizer_cli_hamronized(database):
    file = f'abricate.{database}.tsv'
    perform_cli_test('abricate', file, 'hamronized', db=database)

@pytest.mark.parametrize('database', ['argannot', 'megares', 'ncbi'])
def test_abricate_normalizer_cli_raw(database):
    file = f'abricate.{database}.tsv'
    perform_cli_test('abricate', file, 'raw', db=database)

@pytest.mark.parametrize('folder', ['hamronized', 'raw'])
@pytest.mark.parametrize('mode', ['reads', 'orfs'])
def test_resfinder_cli(folder, mode):
    file = f'resfinder.resfinder.{mode}.tsv'
    perform_cli_test('resfinder', file, folder)

@pytest.mark.parametrize('folder', ['hamronized', 'raw'])
def test_amrfinderplus_cli(folder):
    file = 'amrfinderplus.ncbi.orfs.tsv'
    perform_cli_test('amrfinderplus', file, folder)
