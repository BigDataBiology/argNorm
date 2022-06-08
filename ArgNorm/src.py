"""Source code for the Normalizer class"""

from lib2to3.pytree import Base
import os
import pandas as pd


_ROOT = os.path.abspath(os.path.dirname(__file__))

def get_data(path):
    return os.path.join(_ROOT, 'data', path)

supported_dbs = ['ncbi', 'resfinder', 'sarg', 'deeparg', 'megares', 'argannot']


class HamronizedNormalizer:
    """Major class to perform the normalization"""

    def __init__(self, db):
        self.database = db
        self._set_mapping_table()
        self._aro_col = 'ARO'
        self._set_gene_name_cols()
        self._set_input_gene_name_process_funcs()
        self._set_AROmap_gene_name_process_funcs()
        self._set_input_loading_funcs()
        self._temp_gene_name_col = 'temp_gene_name'

    def _set_mapping_table(self):
        suffix = '_ARO_mapping.tsv'
        self.mapping_tables = dict(
            ncbi=get_data('ncbi' + suffix),
            deeparg=get_data('deeparg' + suffix),
            resfinder=get_data('resfinder' + suffix),
            sarg=get_data('sarg' + suffix),
            megares=get_data('megares' + suffix),
            argannot=get_data('argannot' + suffix)
        )

    def _set_gene_name_cols(self):
        self._gene_name_cols = dict(
            ncbi='gene_symbol',
            deeparg='gene_name',
            resfinder='gene_symbol',
            sarg='gene_symbol',
            megares='reference_accession',
            argannot='reference_accession'
        )
        
    def _set_input_gene_name_process_funcs(self):
        self._preprocess_input_gene_name_funcs = dict(
            ncbi=lambda x: x,
            deeparg=lambda x: x,
            resfinder=lambda x: x,
            sarg=lambda x: x,
            megares=lambda x: x,
            argannot=lambda x: x
        )

    def _set_AROmap_gene_name_process_funcs(self):
        self._preprocess_AROmap_gene_name_funcs = dict(
            ncbi=lambda x: x.split('|')[5],
            deeparg=lambda x: x,
            resfinder=lambda x: '_'.join(x.split('_')[:-1]),
            sarg=lambda x: x,
            megares=lambda x: x.split('|')[0],
            argannot=lambda x: x.split('~~~')[-1]
        )

    def _set_input_loading_funcs(self):
        fun = lambda file: pd.read_csv(file, sep='\t', index_col=0)
        self._input_loading_funcs = dict(
            ncbi=fun,
            deeparg=fun,
            resfinder=fun,
            sarg=fun,
            megares=fun,
            argannot=fun
        )

    def run(self, input_file):
        """Performs the normalization"""
        original_annot = self._load_input(input_file)
        original_annot = self._preprocess_input(original_annot)
        ##
        mapping_table = self.load_mapping_table()
        names = self._preprocess_gene_names(mapping_table['Original ID']).str.lower()
        aros = 'ARO:' + mapping_table[self._aro_col].astype(str).apply(lambda x: x.split('.')[0])
        hashmap = self._make_hashmap(names, aros)
        ##
        original_annot[self._aro_col] = original_annot[self._temp_gene_name_col].\
            str.lower().map(hashmap)
        return original_annot.drop(columns=[self._temp_gene_name_col])

    def _preprocess_input(self, input_annot):
        gene_name_col = self._gene_name_cols[self.database]
        single_gene_name_process = self._preprocess_input_gene_name_funcs[self.database]
        input_annot[self._temp_gene_name_col] = input_annot[gene_name_col].\
            apply(single_gene_name_process)
        return input_annot

    def _preprocess_gene_names(self, gene_names):
        return gene_names.apply(self._preprocess_AROmap_gene_name_funcs[self.database])

    def load_mapping_table(self):
        file = self.mapping_tables[self.database]
        return pd.read_csv(file, sep='\t', index_col=0)

    def _load_input(self, file):
        return self._input_loading_funcs[self.database](file)

    def _make_hashmap(self, keys, vals):
        return dict(zip(keys, vals))


class RawNormalizer(HamronizedNormalizer):

    def __init__(self, db):
        HamronizedNormalizer.__init__(self, db=db)

    def _set_gene_name_cols(self):
        self._gene_name_cols = dict(
            ncbi='GENE',
            deeparg='best-hit',
            resfinder='GENE',
            megares='ACCESSION',
            argannot='ACCESSION'
        )

    # def _set_input_gene_name_process_funcs(self):
        # self._preprocess_input_gene_name_funcs = dict(
        #     ncbi=lambda x: x
        # )

    # def _set_AROmap_gene_name_process_funcs(self):
        # self._preprocess_AROmap_gene_name_funcs = dict(
        #     ncbi=lambda x: x.split('|')[5],
        # )

    def _set_input_loading_funcs(self):
        HamronizedNormalizer._set_input_loading_funcs(self)
        self._input_loading_funcs['sarg'] = lambda x: x # TODO to be updated

    def run(self, input_file):
        return HamronizedNormalizer.run(self, input_file=input_file)
