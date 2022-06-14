"""Major script to run in command line"""

import argparse
from .normalizers import ARGSOAPNormalizer, DeepARGNormalizer, AbricateNormalizer


def main():
    """
    Major function to run when running `argnorm` in shell.
    """
    parser = argparse.ArgumentParser(
        description=('The program is designed for normalizing ARG annotation result '
                      'from different ARG annotation tools and databases to resolve '
                      'their differences in gene naming etc.'),
    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('tool', type=str,
                        choices=['argsoap', 'abricate', 'deeparg'],
                        help='The tool you used to do ARG annotation.')
    parser.add_argument('--db', type=str,
                        choices=['sarg', 'ncbi', 'resfinder', 'deeparg', 'megares', 'argannot'],
                        help='The database you used to do ARG annotation.')
    parser.add_argument('--mode', type=str,
                        choices=['reads', 'orfs', 'both'],
                        help='The tool you used to do ARG annotation.')
    parser.add_argument('--hamronized', action='store_true', help='Use this if the input is hamronized (not hamronized by hAMRonization)')
    parser.add_argument('-i', '--input', type=str, help='The annotation result you have.')
    parser.add_argument('-o', '--output', type=str, help='The file to save normalization results.')
    args = parser.parse_args()

    if args.tool == 'argsoap':
        norm = ARGSOAPNormalizer(database=args.db, is_hamronized=args.hamronized, mode=args.mode)
    elif args.tool == 'deeparg':
        norm = DeepARGNormalizer(database=args.db, is_hamronized=args.hamronized, mode=args.mode)
    elif args.tool == 'abricate':
        norm = AbricateNormalizer(database=args.db, is_hamronized=args.hamronized, mode=args.mode)
    else:
        raise ValueError('Please specify a correct tool name.')
    result = norm.run(input_file=args.input)
    print(result)
    prop_unmapped = ((result.ARO == 'ARO:nan').sum() + result.ARO.isna().sum()) / result.shape[0]
    print(f'{round(1 - prop_unmapped, 3):.2%} args mapped.')
    result.to_csv(args.output, sep='\t')
