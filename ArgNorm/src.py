import pandas as pd
import os

_ROOT = os.path.abspath(os.path.dirname(__file__))

def get_data(path):
    return os.path.join(_ROOT, 'data', path)

supported_dbs = ['ncbi', 'resfinder', 'sarg', 'deeparg', 'megares', 'argannot']

class Normalizer(object):

    def __init__(self, db):
        self.database = db
        suffix = '_ARO_mapping.tsv'
        self.mapping_tables = dict(
            ncbi=get_data('ncbi' + suffix),
            deeparg=get_data('deeparg' + suffix),
            resfinder=get_data('resfinder' + suffix),
            sarg=get_data('sarg' + suffix),
            megares=get_data('megares' + suffix),
            argannot=get_data('argannot' + suffix)
        )
        self.__aro_col = 'ARO'
        self.__gene_name_cols = dict(
            ncbi='gene_symbol',
            deeparg='gene_name',
            resfinder='gene_symbol',
            sarg='gene_symbol',
            megares='gene_name',
            argannot='gene_symbol'
        )
        self.__preprocess_input_gene_name_funcs = dict(
            ncbi=lambda x: x,
            deeparg=lambda x: x,
            resfinder=lambda x: x,
            sarg=lambda x: x,
            megares=lambda x: x,
            argannot=lambda x: x
        )
        self.__preprocess_AROmap_gene_name_funcs = dict(
            ncbi=lambda x: x.split('|')[5],
            deeparg=lambda x: x,
            resfinder=lambda x: '_'.join(x.split('_')[:-1]),
            sarg=lambda x: x,
            megares=lambda x: ':'.join(x.split('|')[1:]),
            argannot=lambda x: x.split('~~~')[1]
        )
        self.__temp_gene_name_col = 'temp_gene_name'
        self.__supported_dbs = supported_dbs

    def run(self, input_file):
        original_annot = self.__load_input(input_file)
        original_annot = self._preprocess_input(original_annot)
        
        mapping_table = self.load_mapping_table()
        names = self._preprocess_gene_names(mapping_table['Original ID']).str.lower()
        aros = 'ARO:' + mapping_table['ARO'].astype(str).apply(lambda x: x.split('.')[0])
        hashmap = self.__make_hashmap(names, aros)
        
        original_annot[self.__aro_col] = original_annot[self.__temp_gene_name_col].str.lower().map(hashmap)
        return original_annot.drop(columns=[self.__temp_gene_name_col])

    def _preprocess_input(self, input_annot):
        gene_name_col = self.__gene_name_cols[self.database]
        single_gene_name_process = self.__preprocess_input_gene_name_funcs[self.database]
        input_annot[self.__temp_gene_name_col] = input_annot[gene_name_col].apply(single_gene_name_process)
        return input_annot

    def _preprocess_gene_names(self, gene_names):
        return gene_names.apply(self.__preprocess_AROmap_gene_name_funcs[self.database])

    def load_mapping_table(self):
        file = self.mapping_tables[self.database]
        return pd.read_csv(file, sep='\t', index_col=0)

    def __load_input(self, file):
        return pd.read_csv(file, sep='\t', index_col=0)

    def __make_hashmap(self, keys, vals):
        return dict(zip(keys, vals))
