import pytest
import argnorm.normalize as argnorm
import pandas as pd
import os

def get_normed(normalizer, input_path):
    normed = normalizer.run(input_path)
    assert not (normed['ARO'] == 'ARO:nan').any()
    return normed


def get_golden_file(category, file_name):
    path = os.path.join('./outputs/', category, file_name)
    return pd.read_csv(path, sep='\t', skiprows=1)

def test_argsoap_normalizer():
    normalizer = argnorm.ARGSOAPNormalizer()
    input_path = './examples/raw/args-oap.sarg.reads.tsv'

    normed = get_normed(normalizer, input_path)
    golden_file = get_golden_file('raw', 'args-oap.sarg.reads.tsv')
    # The numbers in the heading of the raw normed df have the type int, while numbers in heading of golden file have type str
    normed.columns = golden_file.columns

    pd.testing.assert_frame_equal(normed, golden_file)

def test_deeparg_normalizer():
    normalizer = argnorm.DeepARGNormalizer()
    input_path = './examples/raw/deeparg.deeparg.orfs.tsv'

    normed = get_normed(normalizer, input_path)
    golden_file = get_golden_file('raw', 'deeparg.deeparg.orfs.tsv')

    pd.testing.assert_frame_equal(normed, golden_file)


@pytest.mark.parametrize('database', ['argannot', 'megares', 'ncbi', 'resfinder'])
def test_abricate_normalizer(database):
    normalizer = argnorm.AbricateNormalizer(database=database)
    input_path = f'./examples/raw/abricate.{database}.tsv'

    normed = get_normed(normalizer, input_path)
    golden_file = get_golden_file('raw', f'abricate.{database}.tsv')

    pd.testing.assert_frame_equal(normed, golden_file)


@pytest.mark.parametrize('database', ['random_db', 'sarg', 'resfinderfg'])
def test_abricate_validation(database):
    with pytest.raises(Exception):
        argnorm.AbricateNormalizer(database=database)


@pytest.mark.parametrize('mode', ['reads', 'orfs'])
def test_resfinder_normalizer(mode):
    normalizer = argnorm.ResFinderNormalizer()
    input_path = f'./examples/raw/resfinder.resfinder.{mode}.tsv'
    
    normed = get_normed(normalizer, input_path)
    golden_file = get_golden_file('raw', f'resfinder.resfinder.{mode}.tsv')

    pd.testing.assert_frame_equal(normed, golden_file)


@pytest.mark.parametrize('version', ['v3.10.30', 'v4'])
def test_amrfinderplus_normalizer(version):
    normalizer = argnorm.AMRFinderPlusNormalizer()
    input_file = f'./examples/raw/amrfinderplus.ncbi.orfs.{version}.tsv'

    normed = get_normed(normalizer, input_file)
    golden_file = get_golden_file('raw', f'amrfinderplus.ncbi.orfs.{version}.tsv')

    pd.testing.assert_frame_equal(normed, golden_file)


@pytest.mark.parametrize('database', ['argannot', 'resfinder', 'card', 'groot-db', 'groot-core-db'])
def test_groot_normalizer(database):
    if 'groot' not in database:
        groot_db = f'groot-{database}'
    else:
        groot_db = database

    normalizer = argnorm.GrootNormalizer(database=groot_db)
    input_path = f'./examples/raw/groot.{database}.tsv'
    normed = get_normed(normalizer, input_path)
    normed.columns = normed.columns.astype(str)
    golden_file = get_golden_file('raw', f'groot.{database}.tsv')
    pd.testing.assert_frame_equal(normed, golden_file)


def test_hamronization_normalizer():
    normalizer = argnorm.HamronizationNormalizer()

    for file in os.listdir('./examples/hamronized/'):
        if 'abricate.card' in file or 'args-oap.sarg.orfs' in file:
            continue
        if file == 'combined_hamronization_full.tsv':
            normalizer = argnorm.HamronizationNormalizer(skip_on_unsupported_tool=True)

        normed = get_normed(normalizer, f'./examples/hamronized/{file}')
        golden_file = get_golden_file('hamronized', file)
        pd.testing.assert_frame_equal(normed, golden_file)

def test_raw_matching_hamronization():
    files = [
        'resfinder.resfinder.orfs.tsv',
        'resfinder.resfinder.reads.tsv',
        'groot.resfinder.tsv',
        'groot.groot-db.tsv',
        'groot.groot-core-db.tsv',
        'groot.card.tsv',
        'groot.argannot.tsv',
        'amrfinderplus.ncbi.orfs.v3.10.30.tsv'
    ]
    
    hamronization_normalizer = argnorm.HamronizationNormalizer(skip_on_unsupported_tool=True)
    
    for file in files:
        hamronized_file = f'./examples/hamronized/{file}'
        raw_file = f'./examples/raw/{file}'
        
        print(raw_file)
        
        if 'resfinder' in file and not 'groot' in file:
            normalizer = argnorm.ResFinderNormalizer()
        elif 'groot' in file:
            normalizer = argnorm.GrootNormalizer(database=((f"groot-{file.split('.')[1]}") if not 'groot' in file.split('.')[1] else file.split('.')[1]))
        elif 'amrfinderplus' in file:
            normalizer = argnorm.AMRFinderPlusNormalizer()
        else:
            raise Exception('Unable to determine normalizer')        
        
        ham_normed = get_normed(hamronization_normalizer, hamronized_file)
        raw_normed = get_normed(normalizer, raw_file)
        
        pd.testing.assert_series_equal(ham_normed['ARO'], raw_normed['ARO'])