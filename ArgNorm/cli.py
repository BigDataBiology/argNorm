"""Major script to run in command line"""

import argparse
from .src import HamronizedNormalizer, RawNormalizer


def main():
    """
    Major function to run when running `argnorm` in shell.
    """
    parser = argparse.ArgumentParser(
        description=('The program is designed for normalizing ARG annotation result '
                      'from different ARG annotation tools and databases to resolve '
                      'their differences in gene naming etc.'),
    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('database', type=str,
                        choices=['ncbi', 'resfinder', 'deeparg', 'megares', 'argannot'],
                        help='The database you used to do ARG annotation.')
    parser.add_argument('--raw', action='store_true', help='Use this if the input is raw (not hamronized by hAMRonization)')
    parser.add_argument('-i', '--input', type=str, help='The annotation result you have.')
    parser.add_argument('-o', '--output', type=str, help='The file to save normalization results.')
    args = parser.parse_args()
    if not args.raw:
        norm = HamronizedNormalizer(db=args.database)
    else:
        norm = RawNormalizer(db=args.database)
    result = norm.run(input_file=args.input)
    prop_unmapped = ((result.ARO == 'ARO:nan').sum() + result.ARO.isna().sum()) / result.shape[0]
    print(f'{round(1 - prop_unmapped, 3):.2%} args mapped.')
    result.to_csv(args.output, sep='\t')
