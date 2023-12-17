import pytest
import argnorm.normalizers as argnorm


@pytest.mark.parametrize('uses_manual_curation', [True, False])
def test_add_aro_column_argsoap_reads_hamronized(uses_manual_curation):
    norm = argnorm.ARGSOAPNormalizer(is_hamronized=True, mode='reads', uses_manual_curation=uses_manual_curation)
    normed = norm.run(input_file=f'examples/hamronized/args-oap.sarg.reads.tsv')

    assert 'ARO' in normed.columns
    assert set(normed.loc[normed['gene_symbol'] == 'ykkD', 'ARO'].tolist()) == {'ARO:3003064'}

    if uses_manual_curation:
        assert set(normed.loc[normed['gene_symbol'] == 'tetR', 'ARO'].tolist()) == {'ARO:gi|504416866|ref|WP_014603968.1|', 'ARO:gi|553755386|ref|WP_023088115.1|'}
    else:
        assert set(normed.loc[normed['gene_symbol'] == 'tetR', 'ARO'].tolist()) == {'ARO:nan'}


@pytest.mark.parametrize('uses_manual_curation', [True, False])
def test_add_aro_column_argsoap_reads_raw(uses_manual_curation):
    norm = argnorm.ARGSOAPNormalizer(is_hamronized=False, mode='reads', uses_manual_curation=uses_manual_curation)
    normed = norm.run(input_file=f'examples/raw/args-oap.sarg.reads.tsv')

    assert 'ARO' in normed.columns
    assert set(normed.loc[normed[1] == 'gi|489421661|ref|WP_003327389.1|', 'ARO'].tolist()) == {'ARO:3003064'}


@pytest.mark.parametrize('uses_manual_curation', [True, False])
def test_add_aro_column_argsoap_orfs_hamronized(uses_manual_curation):
    norm = argnorm.ARGSOAPNormalizer(is_hamronized=True, mode='orfs', uses_manual_curation=uses_manual_curation)
    normed = norm.run(input_file=f'examples/hamronized/args-oap.sarg.orfs.tsv')
    assert 'ARO' in normed.columns
    assert set(normed.loc[normed['gene_symbol'] == 'ykkD', 'ARO'].tolist()) == {"{'ARO:3003064'}"}

    if uses_manual_curation:
        assert set(normed.loc[normed['gene_symbol'] == 'tetR', 'ARO'].tolist()) == {"{'ARO:0000051'}"}
        assert set(normed.loc[normed['gene_symbol'] == 'mexT', 'ARO'].tolist()) == {"{'ARO:3000814'}"}
    else:
        assert set(normed.loc[normed['gene_symbol'] == 'tetR', 'ARO'].tolist()) == {"{'ARO:3001805', 'ARO:3000559', 'ARO:3003710', 'ARO:nan'}"}
        assert set(normed.loc[normed['gene_symbol'] == 'mexT', 'ARO'].tolist()) == {"{'ARO:nan'}"}


@pytest.mark.parametrize('uses_manual_curation', [True, False])
def test_add_aro_column_argsoap_orfs_raw(uses_manual_curation):
    norm = argnorm.ARGSOAPNormalizer(is_hamronized=False, mode='orfs', uses_manual_curation=uses_manual_curation)
    normed = norm.run(input_file=f'examples/raw/args-oap.sarg.orfs.tsv')
    assert 'ARO' in normed.columns
    assert set(normed.loc[normed[0].str.contains('multidrug__ykkD_train_msa'), 'ARO'].tolist()) == {"{'ARO:3003064'}"}

    if uses_manual_curation:
        assert set(normed.loc[normed[0].str.contains('tetracycline__tetR_train_msa'), 'ARO'].tolist()) == {"{'ARO:0000051'}"}
        assert set(normed.loc[normed[0].str.contains('multidrug__mexT_train_msa'), 'ARO'].tolist()) == {"{'ARO:3000814'}"}
    else:
        assert set(normed.loc[normed[0].str.contains('tetracycline__tetR_train_msa'), 'ARO'].tolist()) == {"{'ARO:3001805', 'ARO:3000559', 'ARO:3003710', 'ARO:nan'}"}
        assert set(normed.loc[normed[0].str.contains('multidrug__mexT_train_msa'), 'ARO'].tolist()) == {"{'ARO:nan'}"}


@pytest.mark.parametrize('uses_manual_curation', [True, False])
def test_add_aro_column_deeparg_hamronized(uses_manual_curation):
    norm = argnorm.DeepARGNormalizer(is_hamronized=True, uses_manual_curation=uses_manual_curation)
    normed = norm.run(input_file=f'examples/hamronized/deeparg.deeparg.orfs.tsv')
    assert 'ARO' in normed.columns
    assert set(normed.loc[normed['gene_symbol'] == 'YKKD', 'ARO'].tolist()) == {'ARO:3003064'}

    if uses_manual_curation:
        assert set(normed.loc[normed['gene_symbol'] == 'MGRB', 'ARO'].tolist()) == {'ARO:mgrB'}
    else:
        assert set(normed.loc[normed['gene_symbol'] == 'MGRB', 'ARO'].tolist()) == {'ARO:nan'}


@pytest.mark.parametrize('uses_manual_curation', [True, False])
def test_add_aro_column_deeparg_raw(uses_manual_curation):
    norm = argnorm.DeepARGNormalizer(is_hamronized=False, uses_manual_curation=uses_manual_curation)
    normed = norm.run(input_file=f'examples/raw/deeparg.deeparg.orfs.tsv')

    assert 'ARO' in normed.columns
    assert normed.set_index('#ARG').loc['YKKD', 'ARO'] == 'ARO:3003064'

    if uses_manual_curation:
        assert list(normed.set_index('#ARG').loc['PENA', 'ARO']) == ['ARO:3004832', 'ARO:3004832']
        assert list(normed.set_index('#ARG').loc['PBP-2X', 'ARO']) == ['ARO:PBP-2X', 'ARO:PBP-2X']
        assert normed.set_index('#ARG').loc['TETR', 'ARO'] == 'ARO:3003479'
        assert normed.set_index('#ARG').loc['PATA', 'ARO'] == 'ARO:3000024'
    else:
        assert list(normed.set_index('#ARG').loc['PENA', 'ARO']) == ['ARO:nan', 'ARO:nan']
        assert list(normed.set_index('#ARG').loc['PBP-2X', 'ARO']) == ['ARO:nan', 'ARO:nan']
        assert normed.set_index('#ARG').loc['TETR', 'ARO'] == 'ARO:nan'
        assert normed.set_index('#ARG').loc['PATA', 'ARO'] == 'ARO:nan'


@pytest.mark.parametrize("database", ['argannot', 'megares', 'ncbi', 'resfinder'])
@pytest.mark.parametrize('uses_manual_curation', [True, False])
def test_add_aro_column_abricate_hamronized(database, uses_manual_curation):
    norm = argnorm.AbricateNormalizer(database=database, is_hamronized=True, uses_manual_curation=uses_manual_curation)
    normed = norm.run(input_file=f'examples/hamronized/abricate.{database}.tsv')

    assert 'ARO' in normed.columns
    assert normed.set_index('input_sequence_id').loc['GMGC10.017_618_532.GPT', 'ARO'] == 'ARO:3001305'

@pytest.mark.parametrize("database", ['argannot', 'megares', 'ncbi'])
@pytest.mark.parametrize('uses_manual_curation', [True, False])
def test_add_aro_column_abricate_raw(database, uses_manual_curation):
    norm = argnorm.AbricateNormalizer(database=database, is_hamronized=False, uses_manual_curation=uses_manual_curation)
    normed = norm.run(input_file=f'examples/raw/abricate.{database}.tsv')

    assert 'ARO' in normed.columns
    assert normed.set_index('SEQUENCE').loc['GMGC10.034_105_239.FOLA', 'ARO'] == 'ARO:3002858'

@pytest.mark.parametrize("mode", ['reads', 'orfs'])
@pytest.mark.parametrize('uses_manual_curation', [True, False])
def test_add_aro_column_resfinder_raw(mode, uses_manual_curation):
    norm = argnorm.ResFinderNormalizer(is_hamronized=False, uses_manual_curation=uses_manual_curation)
    normed = norm.run(input_file=f'examples/raw/resfinder.resfinder.{mode}.tsv')
    assert 'ARO' in normed.columns
    assert normed.set_index('Resistance gene').loc["aph(3')-III", 'ARO'] == 'ARO:3002647'

# Will be getting hamronized inputs from PR#11
@pytest.mark.parametrize('uses_manual_curation', [True, False])
def test_add_aro_column_amrfinder_plus_raw(uses_manual_curation):
    norm = argnorm.AMRFinderPlusNormalizer(is_hamronized=False, uses_manual_curation=uses_manual_curation)
    normed = norm.run(input_file=f'examples/raw/amrfinderplus.ncbi.orfs.tsv')

    assert 'ARO' in normed.columns
    assert normed.set_index('Gene symbol').loc["aph(3')-IIIa", 'ARO'] == 'ARO:3002647'

@pytest.mark.parametrize('uses_manual_curation', [True, False])
def test_add_aro_column_amrfinder_plus_hamronized(uses_manual_curation):
    norm = argnorm.AMRFinderPlusNormalizer(is_hamronized=True, uses_manual_curation=uses_manual_curation)
    normed = norm.run(input_file=f'examples/hamronized/amrfinderplus.ncbi.orfs.tsv')

    assert 'ARO' in normed.columns
    assert normed.set_index('gene_symbol').loc["aph(3')-IIIa", 'ARO'] == 'ARO:3002647'