import pytest
import argNorm.normalizers as argnorm

@pytest.mark.parametrize("is_hamronized", [True, False])
@pytest.mark.parametrize("mode", ['reads', 'orfs'])
def test_add_aro_column_argsoap(is_hamronized, mode):
    example_dir = 'hamronized' if is_hamronized else 'raw'
    norm = argnorm.ARGSOAPNormalizer(is_hamronized=is_hamronized, mode=mode)
    normed = norm.run(input_file=f'examples/{example_dir}/args-oap.sarg.{mode}.tsv')
    assert 'ARO' in normed.columns


@pytest.mark.parametrize("is_hamronized", [True, False])
def test_add_aro_column_deeparg(is_hamronized):
    example_dir = 'hamronized' if is_hamronized else 'raw'
    norm = argnorm.DeepARGNormalizer(is_hamronized=is_hamronized)
    normed = norm.run(input_file=f'examples/{example_dir}/deeparg.deeparg.orfs.tsv')
    assert 'ARO' in normed.columns


@pytest.mark.parametrize("is_hamronized", [True, False])
@pytest.mark.parametrize("db", ['argannot', 'megares', 'ncbi', 'resfinder'])
def test_add_aro_column_abricate(is_hamronized, db):
    if not is_hamronized and db == 'resfinder':
        pass  
        # TODO We are missing this example test input.
    else:
        example_dir = 'hamronized' if is_hamronized else 'raw'
        norm = argnorm.AbricateNormalizer(database=db, is_hamronized=is_hamronized)
        normed = norm.run(input_file=f'examples/{example_dir}/abricate.{db}.tsv')
        assert 'ARO' in normed.columns
