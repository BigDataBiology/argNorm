import pytest
import argnorm.normalizers as argnorm

def test_argsoap_normalizer_hamronized():
    norm = argnorm.ARGSOAPNormalizer(is_hamronized=True)
    normed = norm.run(input_file=f'examples/hamronized/args-oap.sarg.reads.tsv')
    
    assert 'ARO' in normed.columns
    assert set(normed.loc[normed['gene_symbol'] == 'ykkD', 'ARO'].tolist()) == {'ARO:3003064'}

def test_argsoap_normalizer_raw():
    norm = argnorm.ARGSOAPNormalizer(is_hamronized=False)
    normed = norm.run(input_file=f'examples/raw/args-oap.sarg.reads.tsv')
    
    assert 'ARO' in normed.columns
    assert set(normed.loc[normed[1] == 'gi|489421661|ref|WP_003327389.1|', 'ARO'].tolist()) == {'ARO:3003064'}
    
def test_deeparg_hamronized():
    norm = argnorm.DeepARGNormalizer(is_hamronized=True)
    normed = norm.run(input_file=f'examples/hamronized/deeparg.deeparg.orfs.tsv')
    
    assert 'ARO' in normed.columns
    assert set(normed.loc[normed['gene_symbol'] == 'YKKD', 'ARO'].tolist()) == {'ARO:3003064'}
    
def test_deeparg_raw():
    norm = argnorm.DeepARGNormalizer(is_hamronized=False)
    normed = norm.run(input_file=f'examples/raw/deeparg.deeparg.orfs.tsv')
    
    assert 'ARO' in normed.columns
    assert normed.set_index('#ARG').loc['YKKD', 'ARO'] == 'ARO:3003064'
    assert list(normed.set_index('#ARG').loc['PENA', 'ARO']) == ['ARO:3003042', 'ARO:3003042']
    assert list(normed.set_index('#ARG').loc['PBP-2X', 'ARO']) == ['ARO:3003937', 'ARO:3003937']
    assert normed.set_index('#ARG').loc['TETR', 'ARO'] == 'ARO:3003479'

@pytest.mark.parametrize("database", ['argannot', 'megares', 'ncbi', 'resfinder'])
def test_abricate_hamronized(database):
    norm = argnorm.AbricateNormalizer(database=database, is_hamronized=True)
    normed = norm.run(input_file=f'examples/hamronized/abricate.{database}.tsv')

    assert 'ARO' in normed.columns
    assert normed.set_index('input_sequence_id').loc['GMGC10.017_618_532.GPT', 'ARO'] == 'ARO:3001305'

@pytest.mark.parametrize("database", ['argannot', 'megares', 'ncbi'])
def test_abricate_raw(database):
    norm = argnorm.AbricateNormalizer(database=database, is_hamronized=False)
    normed = norm.run(input_file=f'examples/raw/abricate.{database}.tsv')

    assert 'ARO' in normed.columns
    assert normed.set_index('SEQUENCE').loc['GMGC10.034_105_239.FOLA', 'ARO'] == 'ARO:3002858'
    
@pytest.mark.parametrize("mode", ["reads", "orfs"])
def test_resfinder_hamronized(mode):
    norm = argnorm.ResFinderNormalizer(is_hamronized=True)
    normed = norm.run(input_file=f'examples/hamronized/resfinder.resfinder.{mode}.tsv')
    
    assert 'ARO' in normed.columns
    assert normed.set_index('gene_symbol').loc["aph(3')-III", 'ARO'] == 'ARO:3002647'

@pytest.mark.parametrize("mode", ['reads', 'orfs'])
def test_resfinder_raw(mode):
    norm = argnorm.ResFinderNormalizer(is_hamronized=False)
    normed = norm.run(input_file=f'examples/raw/resfinder.resfinder.{mode}.tsv')
    
    assert 'ARO' in normed.columns
    assert normed.set_index('Resistance gene').loc["aph(3')-III", 'ARO'] == 'ARO:3002647'
    
def test_amrfinderplus_hamronized():
    norm = argnorm.AMRFinderPlusNormalizer(is_hamronized=True)
    normed = norm.run(input_file=f'examples/hamronized/amrfinderplus.ncbi.orfs.tsv')
    
    assert 'ARO' in normed.columns
    assert normed.set_index('gene_symbol').loc["aph(3')-IIIa", 'ARO'] == 'ARO:3002647'
    
def test_amrfinderplus_raw():
    norm = argnorm.AMRFinderPlusNormalizer(is_hamronized=False)
    normed = norm.run(input_file=f'examples/raw/amrfinderplus.ncbi.orfs.tsv')
    
    assert 'ARO' in normed.columns
    assert normed.set_index('Gene symbol').loc["aph(3')-IIIa", 'ARO'] == "ARO:3002647"
