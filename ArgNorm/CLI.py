import sys
import os
import argparse
from . import Normalizer


def main():
    parser = argparse.ArgumentParser(
        description=('The program is designed to help you to normalize the annotation result \
                      by different ARG annotation tools and databases to resolve their differences in gene naming etc.'),
		formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('database', type=str, choices=['ncbi', 'resfinder', 'sarg', 'deeparg'],
        help='The database you used to do ARG annotation.')
    parser.add_argument('input', type=str, help='The annotation result you have.')
    parser.add_argument('output', type=str, help='The file to save normalization results.')
    