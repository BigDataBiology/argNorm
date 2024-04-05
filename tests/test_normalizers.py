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
 
@pytest.mark.parametrize('hamronized', [True, False])
def test_argsoap_normalizer(hamronized):
    folder = 'hamronized' if hamronized else 'raw'
    normalizer = argnorm.ARGSOAPNormalizer(is_hamronized=hamronized)
    input_path = f'./examples/{folder}/args-oap.sarg.reads.tsv'

    normed = get_normed(normalizer, input_path)
    golden_file = pd.read_csv(os.path.join('./outputs/', folder, 'args-oap.sarg.reads.tsv'), sep='\t')
    # The numbers in the heading of the raw normed df have the type int, while numbers in heading of golden file have type str
    normed.columns = golden_file.columns

    pd.testing.assert_frame_equal(normed, golden_file)

@pytest.mark.parametrize('hamronized', [True, False])
def test_deeparg_normalizer(hamronized):
    folder = 'hamronized' if hamronized else 'raw'
    normalizer = argnorm.DeepARGNormalizer(is_hamronized=hamronized)
    input_path = f'./examples/{folder}/deeparg.deeparg.orfs.tsv'

    normed = get_normed(normalizer, input_path)
    golden_file = pd.read_csv(os.path.join('./outputs/', folder, 'deeparg.deeparg.orfs.tsv'), sep='\t')

    pd.testing.assert_frame_equal(normed, golden_file)


@pytest.mark.parametrize('database', ['argannot', 'megares', 'ncbi', 'resfinder', 'resfinderfg'])
def test_abricate_normalizer_hamronized(database):
    normalizer = argnorm.AbricateNormalizer(database=database, is_hamronized=True)
    input_path = f'./examples/hamronized/abricate.{database}.tsv'

    normed = get_normed(normalizer, input_path)
    golden_file = pd.read_csv(os.path.join('./outputs/', 'hamronized', f'abricate.{database}.tsv'), sep='\t')

    pd.testing.assert_frame_equal(normed, golden_file)

@pytest.mark.parametrize('database', ['argannot', 'megares', 'ncbi'])
def test_abricate_normalizer_raw(database):
    normalizer = argnorm.AbricateNormalizer(database=database, is_hamronized=False)
    input_path = f'./examples/raw/abricate.{database}.tsv'

    normed = get_normed(normalizer, input_path)
    golden_file = pd.read_csv(os.path.join('./outputs/', 'raw', f'abricate.{database}.tsv'), sep='\t')

    pd.testing.assert_frame_equal(normed, golden_file)

def test_abricate_validation_hamronized():
    with pytest.raises(Exception):
        normalizer = argnorm.AbricateNormalizer(database='random_db', is_hamronized=True)

@pytest.mark.parametrize('database', ['random_db', 'sarg', 'resfinderfg'])
def test_abricate_validation_raw(database):
    with pytest.raises(Exception):
        normalizer = argnorm.AbricateNormalizer(database=database, is_hamronized=False)

@pytest.mark.parametrize('hamronized', [True, False])
@pytest.mark.parametrize('mode', ['reads', 'orfs'])
def test_resfinder_normalizer(hamronized, mode):
    folder = 'hamronized' if hamronized else 'raw'
    normalizer = argnorm.ResFinderNormalizer(is_hamronized=hamronized)
    input_path = f'./examples/{folder}/resfinder.resfinder.{mode}.tsv'
    
    normed = get_normed(normalizer, input_path)
    golden_file = pd.read_csv(os.path.join('./outputs/', folder, f'resfinder.resfinder.{mode}.tsv'), sep='\t')

    pd.testing.assert_frame_equal(normed, golden_file)

@pytest.mark.parametrize('hamronized', [True, False])
def test_amrfinderplus_normalizer(hamronized):
    folder = 'hamronized' if hamronized else 'raw'
    normalizer = argnorm.AMRFinderPlusNormalizer(is_hamronized=hamronized)
    input_file = f'./examples/{folder}/amrfinderplus.ncbi.orfs.tsv'
    
    normed = get_normed(normalizer, input_file)
    golden_file = pd.read_csv(os.path.join('./outputs/', folder, f'amrfinderplus.ncbi.orfs.tsv'), sep='\t')

    pd.testing.assert_frame_equal(normed, golden_file)