import pytest
import argNorm.normalizers as argnorm


@pytest.mark.parametrize("is_hamronized", [True, False])
@pytest.mark.parametrize("mode", ['reads', 'orfs'])
def test_run(is_hamronized, mode):
    example_dir = 'hamronized' if is_hamronized else 'raw'

    norm = argnorm.ARGSOAPNormalizer(is_hamronized=is_hamronized, mode=mode)
    norm.run(input_file=f'examples/{example_dir}/args-oap.sarg.{mode}.tsv')
