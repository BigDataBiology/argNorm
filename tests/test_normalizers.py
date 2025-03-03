import pytest
import argnorm.normalize as argnorm
import pandas as pd
import os
import numpy as np
from warnings import warn
    
def get_normed(normalizer, input_path):
    normed = normalizer.run(input_path)
    normed = normed.apply(lambda x: x.str.strip() if isinstance(x, str) else x).replace('', np.nan)
    return normed

def get_golden_file(path):
    return pd.read_csv(path, sep='\t', skiprows=1)

def test_argsoap_normalizer():
    normalizer = argnorm.ARGSOAPNormalizer()
    input_path = './examples/raw/args-oap.sarg.reads.tsv'

    normed = get_normed(normalizer, input_path)
    golden_file = get_golden_file(os.path.join('./outputs/', 'raw', 'args-oap.sarg.reads.tsv'))
    # The numbers in the heading of the raw normed df have the type int, while numbers in heading of golden file have type str
    normed.columns = golden_file.columns

    pd.testing.assert_frame_equal(normed, golden_file)

def test_deeparg_normalizer():
    normalizer = argnorm.DeepARGNormalizer()
    input_path = './examples/raw/deeparg.deeparg.orfs.tsv'

    normed = get_normed(normalizer, input_path)
    golden_file = get_golden_file(os.path.join('./outputs/', 'raw', 'deeparg.deeparg.orfs.tsv'))

    pd.testing.assert_frame_equal(normed, golden_file)

@pytest.mark.parametrize('database', ['argannot', 'megares', 'ncbi'])
def test_abricate_normalizer(database):
    normalizer = argnorm.AbricateNormalizer(database=database)
    input_path = f'./examples/raw/abricate.{database}.tsv'

    normed = get_normed(normalizer, input_path)
    golden_file = get_golden_file(os.path.join('./outputs/', 'raw', f'abricate.{database}.tsv'))

    pd.testing.assert_frame_equal(normed, golden_file)

@pytest.mark.parametrize('database', ['random_db', 'sarg', 'resfinderfg'])
def test_abricate_validation(database):
    with pytest.raises(Exception):
        normalizer = argnorm.AbricateNormalizer(database=database)

@pytest.mark.parametrize('mode', ['reads', 'orfs'])
def test_resfinder_normalizer(mode):
    normalizer = argnorm.ResFinderNormalizer()
    input_path = f'./examples/raw/resfinder.resfinder.{mode}.tsv'
    
    normed = get_normed(normalizer, input_path)
    golden_file = get_golden_file(os.path.join('./outputs/', 'raw', f'resfinder.resfinder.{mode}.tsv'))

    pd.testing.assert_frame_equal(normed, golden_file)

@pytest.mark.parametrize('version', ['v3.10.30', 'v4'])
def test_amrfinderplus_normalizer(version):
    normalizer = argnorm.AMRFinderPlusNormalizer()
    input_file = f'./examples/raw/amrfinderplus.ncbi.orfs.{version}.tsv'
    
    normed = get_normed(normalizer, input_file)
    golden_file = get_golden_file(os.path.join('./outputs/', 'raw', f'amrfinderplus.ncbi.orfs.{version}.tsv'))

    pd.testing.assert_frame_equal(normed, golden_file)

@pytest.mark.parametrize('database', ['argannot', 'resfinder', 'card', 'groot-db', 'groot-core-db'])
def test_groot_normalizer(database):
    if 'groot' not in database:
        normalizer = argnorm.GrootNormalizer(database=f'groot-{database}')
    else:
        normalizer = argnorm.GrootNormalizer(database=f'groot-db')
        
    input_path = f'./examples/raw/groot.{database}.tsv'
    normed = get_normed(normalizer, input_path)
    normed.columns = normed.columns.astype(str)
    golden_file = get_golden_file(os.path.join('./outputs/', 'raw', f'groot.{database}.tsv'))
    pd.testing.assert_frame_equal(normed, golden_file)

def test_hamronization_normalizer():
    normalizer = argnorm.HamronizationNormalizer()

    for file in os.listdir('./examples/hamronized/'):
        if 'abricate.card' in file or 'args-oap.sarg.orfs' in file:
            continue
        if file == 'combined_hamronization_full.tsv':
            normalizer = argnorm.HamronizationNormalizer(skip_on_unsupported_tool=True)

        normed = get_normed(normalizer, f'./examples/hamronized/{file}')
        golden_file = get_golden_file(os.path.join('./outputs', 'hamronized', file))

        pd.testing.assert_frame_equal(normed, golden_file)
