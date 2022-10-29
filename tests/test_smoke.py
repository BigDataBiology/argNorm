import pytest
import argnorm.normalizers as argnorm

@pytest.mark.parametrize("is_hamronized", [True, False])
@pytest.mark.parametrize("mode", ['reads', 'orfs'])
def test_add_aro_column_argsoap(is_hamronized, mode):
    example_dir = 'hamronized' if is_hamronized else 'raw'
    norm = argnorm.ARGSOAPNormalizer(is_hamronized=is_hamronized, mode=mode)
    normed = norm.run(input_file=f'examples/{example_dir}/args-oap.sarg.{mode}.tsv')
    assert 'ARO' in normed.columns
    if is_hamronized and mode == 'orfs':  # argsoap orfs mode predict classes that can be mapped to multiple AROs (must be a set).
        assert set(normed.loc[normed['gene_symbol'] == 'ykkD', 'ARO'].tolist()) == {"{'ARO:3003064'}"}
    elif not is_hamronized and mode == 'orfs':
        assert set(normed.loc[normed[0].str.contains('multidrug__ykkD_train_msa'), 'ARO'].tolist()) == {"{'ARO:3003064'}"}
    elif is_hamronized and mode == 'reads':
        assert set(normed.loc[normed['gene_symbol'] == 'ykkD', 'ARO'].tolist()) == {'ARO:3003064'}
    else:
        assert set(normed.loc[normed[1] == 'gi|489421661|ref|WP_003327389.1|', 'ARO'].tolist()) == {'ARO:3003064'}


@pytest.mark.parametrize("is_hamronized", [True, False])
def test_add_aro_column_deeparg(is_hamronized):
    example_dir = 'hamronized' if is_hamronized else 'raw'
    norm = argnorm.DeepARGNormalizer(is_hamronized=is_hamronized)
    normed = norm.run(input_file=f'examples/{example_dir}/deeparg.deeparg.orfs.tsv')
    assert 'ARO' in normed.columns
    if is_hamronized:
        assert set(normed.loc[normed['gene_symbol'] == 'YKKD', 'ARO'].tolist()) == {'ARO:3003064'}
    else:
        assert normed.set_index('#ARG').loc['YKKD', 'ARO'] == 'ARO:3003064'


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
        if is_hamronized: 
            assert normed.set_index('input_sequence_id' ).loc['GMGC10.017_618_532.GPT', 'ARO'] == 'ARO:3001305'
        else:
            assert normed.set_index('SEQUENCE').loc['GMGC10.034_105_239.FOLA', 'ARO'] == 'ARO:3002858'


@pytest.mark.parametrize("is_hamronized", [False])
@pytest.mark.parametrize("mode", ['reads', 'orfs'])
def test_add_aro_column_resfinder(is_hamronized, mode):  # We only have raw inputs now.
    example_dir = 'raw'
    norm = argnorm.ResFinderNormalizer(is_hamronized=is_hamronized)
    normed = norm.run(input_file=f'examples/{example_dir}/resfinder.resfinder.{mode}.tsv')
    assert 'ARO' in normed.columns
    assert normed.set_index('Resistance gene').loc["aph(3')-III", 'ARO'] == 'ARO:3002647'


@pytest.mark.parametrize("is_hamronized", [False])
@pytest.mark.parametrize("mode", ['orfs'])
def test_add_aro_column_amrfinderplus(is_hamronized, mode):  # We only have raw inputs now.
    example_dir = 'raw'
    norm = argnorm.AMRFinderPlusNormalizer(is_hamronized=is_hamronized)
    normed = norm.run(input_file=f'examples/{example_dir}/amrfinderplus.ncbi.{mode}.tsv')
    assert 'ARO' in normed.columns
    assert normed.set_index('Gene symbol').loc["aph(3')-IIIa", 'ARO'] == 'ARO:3002647'
