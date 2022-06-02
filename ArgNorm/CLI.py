import argparse
from .src import Normalizer
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description=('The program is designed for normalizing ARG annotation result '
                     'from different ARG annotation tools and databases to resolve their differences in gene naming etc.'),
		formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('database', type=str, 
        choices=['ncbi', 'resfinder', 'deeparg', 'megares', 'argannot'],
        help='The database you used to do ARG annotation.')
    parser.add_argument('-i', '--input', type=str, help='The annotation result you have.')
    parser.add_argument('-o', '--output', type=str, help='The file to save normalization results.')
    
    args = parser.parse_args()
    norm = Normalizer(db=args.database)
    result = norm.run(input_file=args.input)
    print('{:.2%} args mapped.'.format(round(1 - ((result.ARO == 'ARO:nan').sum() + result.ARO.isna().sum()) / result.shape[0], 3)))
    result.to_csv(args.output, sep='\t')
