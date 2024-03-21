import pytest
import argnorm.normalize as argnorm
import pandas as pd
import os
import numpy as np

def assert_params(normed, golden_file):
    assert normed['ARO'].equals(golden_file['ARO'])
    assert normed['confers_resistance_to'].equals(golden_file['confers_resistance_to'])
    assert normed['resistance_to_drug_classes'].equals(golden_file['resistance_to_drug_classes'])

def get_normed(normalizer, input_path):
    normed = normalizer.run(input_path)
    normed = normed.apply(lambda x: x.str.strip() if isinstance(x, str) else x).replace('', np.nan)
    return normed

@pytest.mark.parametrize('hamronized', [True, False])
def test_argsoap_normalizer(hamronized):
    folder = 'hamronized' if hamronized else 'raw'
    normalizer = argnorm.ARGSOAPNormalizer(is_hamronized=hamronized)
    input_path = f'./testing/examples/{folder}/args-oap.sarg.reads.tsv'

    normed = get_normed(normalizer, input_path)
    golden_file = pd.read_csv(os.path.join('./testing/outputs/', folder, 'args-oap.sarg.reads.tsv'), sep='\t')

    assert_params(normed, golden_file)

@pytest.mark.parametrize('hamronized', [True, False])
def test_deeparg_normalizer(hamronized):
    folder = 'hamronized' if hamronized else 'raw'
    normalizer = argnorm.DeepARGNormalizer(is_hamronized=hamronized)
    input_path = f'./testing/examples/{folder}/deeparg.deeparg.orfs.tsv'

    normed = get_normed(normalizer, input_path)
    golden_file = pd.read_csv(os.path.join('./testing/outputs/', folder, 'deeparg.deeparg.orfs.tsv'), sep='\t')

    assert_params(normed, golden_file)

@pytest.mark.parametrize('database', ['argannot', 'megares', 'ncbi', 'resfinder'])
def test_abricate_normalizer_hamronized(database):
    normalizer = argnorm.AbricateNormalizer(database=database, is_hamronized=True)
    input_path = f'./testing/examples/hamronized/abricate.{database}.tsv'

    normed = get_normed(normalizer, input_path)
    golden_file = pd.read_csv(os.path.join('./testing/outputs/', 'hamronized', f'abricate.{database}.tsv'), sep='\t')

    assert_params(normed, golden_file)

@pytest.mark.parametrize('database', ['argannot', 'megares', 'ncbi'])
def test_abricate_normalizer_raw(database):
    normalizer = argnorm.AbricateNormalizer(database=database, is_hamronized=False)
    input_path = f'./testing/examples/raw/abricate.{database}.tsv'
    
    normed = get_normed(normalizer, input_path)
    golden_file = pd.read_csv(os.path.join('./testing/outputs/', 'raw', f'abricate.{database}.tsv'), sep='\t')

    assert_params(normed, golden_file)

@pytest.mark.parametrize('hamronized', [True, False])
@pytest.mark.parametrize('mode', ['reads', 'orfs'])
def test_resfinder_normalizer(hamronized, mode):
    folder = 'hamronized' if hamronized else 'raw'
    normalizer = argnorm.ResFinderNormalizer(is_hamronized=hamronized)
    input_path = f'./testing/examples/{folder}/resfinder.resfinder.{mode}.tsv'
    
    normed = get_normed(normalizer, input_path)
    golden_file = pd.read_csv(os.path.join('./testing/outputs/', folder, f'resfinder.resfinder.{mode}.tsv'), sep='\t')

    assert_params(normed, golden_file)

@pytest.mark.parametrize('hamronized', [True, False])
def test_amrfinderplus_normalizer(hamronized):
    folder = 'hamronized' if hamronized else 'raw'
    normalizer = argnorm.AMRFinderPlusNormalizer(is_hamronized=hamronized)
    input_file = f'./testing/examples/{folder}/amrfinderplus.ncbi.orfs.tsv'
    
    normed = get_normed(normalizer, input_file)
    golden_file = pd.read_csv(os.path.join('./testing/outputs/', folder, f'amrfinderplus.ncbi.orfs.tsv'), sep='\t')

    assert_params(normed, golden_file)