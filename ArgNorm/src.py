from . import get_data
import pandas as pd


class Normalizer(object):

    def __init__(self, db):
        self.database = db
        suffix = '_ARO_mapping.tsv'
        self.mapping_tables = dict(
            ncbi=get_data('data/ncbi' + suffix),
            deeparg=get_data('data/deeparg' + suffix),
            resfinder=get_data('data/resfinder' + suffix),
            sarg=get_data('data/sarg' + suffix)
        )
        self.aro_col = 'ARO'
        self.gene_name_col = 'Original ID'

    def preprocess(self):
        pass

    def load_mapping_table(self):
        file = self.mapping_tables[self.database]
        return pd.read_csv(file, sep='\t', index_col=0)

    def __load_input(self, file):
        return pd.read_csv(file, sep='\t', index_col=0)

    def __make_hashmap(self, keys, vals):
        return dict(zip(keys, vals))

    def run(self, input_file, output_file):
        input = self.__load_input(input_file)
        mapping_table = self.load_mapping_table()
        names = mapping_table[self.gene_name_col]
        aros = mapping_table[self.aro_col]
        hashmap = self.__make_hashmap(names, aros)
        result = input.map(hashmap)
        # save...