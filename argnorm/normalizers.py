import pandas as pd
# Pandas <1.5 does not have copy_on_write option
try:
    pd.options.mode.copy_on_write = True
except pd.errors.OptionError:
    pass
from .drug_categorization import confers_resistance_to, drugs_to_drug_classes
from .lib import get_aro_mapping_table, get_aro_ontology
from .lib import MAPPING_TABLE_ARO_COL, TARGET_ARO_COL, DATABASES, CUT_OFF_COL
import sys
from warnings import warn


def _confers_resistance_to_or_nan(aro_id):
    import numpy as np
    r = confers_resistance_to(aro_id)
    if not r:
        return np.nan
    return r


class BaseNormalizer:
    """
    Inherit this class and customize subclass methods to implement the normalization of new databases/formats.
    """

    def __init__(self, database=None) -> None:
        self.database = database

    def run(self, input_file : str):
        """
        Main normalization pipeline.
        """
        original_annot = self.load_input(input_file)
        input_genes = self.get_input_ids(original_annot)

        aro_table = self.load_mapping_table()
        aro_table.set_index(self.preprocess_ref_genes(
            aro_table.index
        ), inplace=True)

        mapping = aro_table[MAPPING_TABLE_ARO_COL].to_dict()
        cut_offs = aro_table[CUT_OFF_COL].to_dict()

        ARO = get_aro_ontology()

        _aro_name_cache = {}
        def get_aro_name(aro_id):
            if aro_id not in _aro_name_cache:
                _aro_name_cache[aro_id] = ARO[aro_id].name
            return _aro_name_cache[aro_id]

        def aro_id_to_names(aro_ids):
            return ','.join([get_aro_name(aro_id) for aro_id in aro_ids])

        original_annot[TARGET_ARO_COL] = input_genes.map(mapping)
        original_annot['ARO_name'] = original_annot[TARGET_ARO_COL].map(
                get_aro_name,
                na_action='ignore')
        original_annot[CUT_OFF_COL] = input_genes.map(cut_offs)

        confers_resistance = original_annot[TARGET_ARO_COL].map(
                _confers_resistance_to_or_nan,
                na_action='ignore')
        resistance_to_drug_classes = confers_resistance.map(
                drugs_to_drug_classes,
                na_action='ignore')
        original_annot['confers_resistance_to'] = confers_resistance.map(
                ','.join,
                na_action='ignore')
        original_annot['confers_resistance_to_names'] = confers_resistance.map(
                aro_id_to_names,
                na_action='ignore')
        original_annot['resistance_to_drug_classes'] = resistance_to_drug_classes.map(
                ','.join,
                na_action='ignore')
        original_annot['resistance_to_drug_classes_names'] = resistance_to_drug_classes.map(
                aro_id_to_names,
                na_action='ignore')

        return original_annot

    def load_mapping_table(self):
        """
        Load the mapping table for the database, defaults to loading a single database mapping table (defined by self.database).
        """
        return get_aro_mapping_table(self.database)


    def preprocess_ref_genes(self, ref_genes):
        """
        Customize this when ref gene and input gene can not exactly match.
        """
        return ref_genes

    def load_input(self, input_file):
        """
        Customize this when it fails to parse the input data.
        """
        return pd.read_csv(input_file, sep='\t')


class ARGSOAPNormalizer(BaseNormalizer):
    def __init__(self, database=None) -> None:
        super().__init__(database='sarg')

    def get_input_ids(self, itable):
        return itable[1]

    def load_input(self, input_file):
        return pd.read_csv(input_file, sep='\t', header=None)


class DeepARGNormalizer(BaseNormalizer):
    def __init__(self, database=None) -> None:
        super().__init__(database='deeparg')

    def get_input_ids(self, itable):
        return itable['best-hit']


class ResFinderNormalizer(BaseNormalizer):
    def __init__(self, database=None) -> None:
        super().__init__(database='resfinder')

    def get_input_ids(self, itable):
        gene_identifier = 'Resistance gene'
        accession = 'Accession no.'
        return pd.Series(itable[gene_identifier] + '_' + itable[accession])

    @staticmethod
    def preprocess_ref_genes(ref_genes):
        split_genes = ref_genes.str.split('_')
        return pd.Series(split_genes.str[0] + '_' + split_genes.str[-1])


class AMRFinderPlusNormalizer(BaseNormalizer):
    def __init__(self, database=None) -> None:
        super().__init__(database='ncbi')

    def get_input_ids(self, itable):
        if 'Accession of closest sequence' in itable.columns and 'Sequence name' in itable.columns:
            accession = 'Accession of closest sequence'
        elif 'Closest reference accession' in itable.columns and 'Element name' in itable.columns:
            accession = 'Closest reference accession'
        else:
            raise NotImplementedError(
                    'Unsupported AMRFinderPlus version detected. '
                    'Supported/tested versions are v3.10.30 or v4.0.19.')

        return pd.Series(itable[accession])
    
    def load_mapping_table(self):
        protein_id_mapping_table = get_aro_mapping_table(self.database)
        refseq_id_mapping_table = get_aro_mapping_table(self.database)
        
        protein_id_index = pd.Series([gene.split('|')[1] for gene in protein_id_mapping_table.index])
        refseq_id_index = pd.Series([gene.split('|')[2] for gene in refseq_id_mapping_table.index])
                
        protein_id_mapping_table.set_index(protein_id_index, inplace=True)
        refseq_id_mapping_table.set_index(refseq_id_index, inplace=True)
        
        return pd.concat([protein_id_mapping_table, refseq_id_mapping_table])


class AbricateNormalizer(BaseNormalizer):
    def __init__(self, database=None) -> None:
        if database not in ['ncbi', 'deeparg', 'resfinder', 'megares', 'argannot']:
            raise Exception(f'{database} is not a supported database.')

        super().__init__(database)

    def get_input_ids(self, itable):
        col = dict(
            ncbi='ACCESSION',
            deeparg='best-hit',
            resfinder=dict(gene_identifier='GENE', accession='ACCESSION'),
            megares='ACCESSION',
            argannot='ACCESSION'
        )[self.database]

        if self.database == 'resfinderfg':
            return itable[col].str.split('|').str[1]
        if self.database == 'resfinder':
            return pd.Series(itable[col['gene_identifier']] + '_' + itable[col['accession']])
        return itable[col]

    @staticmethod
    def preprocess_argannot_ref_gene(ref_gene):
        split_str = ref_gene.split(':')
        if not str(split_str[2][0]).isnumeric() and '-' not in split_str[2]:
            return ':'.join([split_str[1], split_str[3]])
        return ':'.join(ref_gene.split(':')[1:3])

    def preprocess_ref_genes(self, ref_genes):
        process_funcs_by_db = dict(
            resfinder=lambda x: '_'.join([x.split('_')[0], x.split('_')[1], x.split('_')[-1]]),
            megares=lambda x: x.split('|')[0],
            argannot=self.preprocess_argannot_ref_gene,
            resfinderfg=lambda x: x.split('|')[1]
        )
            
        if self.database in ['database', 'sarg', 'ncbi']:
            return ref_genes
        return ref_genes.map(process_funcs_by_db[self.database])
    
    def load_mapping_table(self):
        if self.database != 'ncbi':
            return super().load_mapping_table()
        else:
            protein_id_mapping_table = get_aro_mapping_table(self.database)
            refseq_id_mapping_table = get_aro_mapping_table(self.database)
            
            protein_id_index = pd.Series([gene.split('|')[1] for gene in protein_id_mapping_table.index])
            refseq_id_index = pd.Series([gene.split('|')[2] for gene in refseq_id_mapping_table.index])
                    
            protein_id_mapping_table.set_index(protein_id_index, inplace=True)
            refseq_id_mapping_table.set_index(refseq_id_index, inplace=True)
            
            return pd.concat([protein_id_mapping_table, refseq_id_mapping_table])


class GrootNormalizer(BaseNormalizer):
    def __init__(self, database=None) -> None:
        if database not in ['groot-argannot', 'groot-resfinder', 'groot-db', 'groot-core-db', 'groot-card']:
            raise NotImplementedError(f'{database} is not a supported database for groot.')
        super().__init__(database)

    def load_input(self, input_file):
        return pd.read_csv(input_file, sep='\t', header=None)

    @staticmethod
    def preprocess_groot_db_inputs(gene_name):
        tokens = str(gene_name).split('__')
        processed_gene_name = tokens[1]
        if 'card' in tokens[0].lower():
            processed_gene_name = processed_gene_name.split('|')[-1]
        return processed_gene_name

    def get_input_ids(self, itable):
        if self.database == 'groot-argannot':
            return itable[0].map(lambda x: x.split('~~~')[-1])

        if self.database == 'groot-card':
            return itable[0].map(lambda x: x.split('.')[0])

        if self.database in ['groot-db', 'groot-core-db']:
            return itable[0].map(self.preprocess_groot_db_inputs)

        return itable[0]

    def preprocess_ref_genes(self, ref_genes):
        if self.database == 'groot-argannot':
            return ref_genes.map(lambda x: ':'.join(str(x).split(':')[1:3]))
        return ref_genes


class HamronizationNormalizer(BaseNormalizer):
    def __init__(self, database=None, skip_on_unsupported_tool=False):
        super().__init__(database)
        self.skip_on_unsupported_tool = skip_on_unsupported_tool

    def get_input_ids(self, original_annot):
        input_id_lookup = {
            'argsoap': lambda x: x['reference_accession'],
            'deeparg': lambda x: x['gene_name'],
            'resfinder': lambda x: x['gene_name'] + '_' + x['reference_accession'],
            'amrfinderplus': lambda x: x['reference_accession'],
            'groot': {
                'groot-card': lambda x: x['gene_name'].split('.')[0],
                'groot-argannot': lambda x: x['gene_name'].split('~~~')[-1],
                'groot-db': GrootNormalizer.preprocess_groot_db_inputs,
                'groot-resfinder': lambda x: x['gene_name']
            },
            'abricate': {
                'sarg': lambda x: x['gene_symbol'],
                'megares': lambda x: x['reference_accession'],
                'argannot': lambda x: x['reference_accession'],
                'resfinderfg': lambda x: x['gene_name'].split('|')[1],
                'deeparg': lambda x: x['gene_name'],
                'resfinder': lambda x: x['gene_name'] + '_' + x['reference_accession'],
                'ncbi': lambda x: x['reference_accession']
            }
        }
        input_genes = []
        for _,row in original_annot.iterrows():
            analysis_software = row['analysis_software_name'].lower() \
                                    .replace('-', '') \
                                    .replace('_', '') \
                                    .replace(' ', '')

            if 'ncbi' in analysis_software:
                analysis_software = 'amrfinderplus'

            for tool in input_id_lookup:
                if tool in analysis_software:
                    analysis_software = tool
                    break
            else:
                if self.skip_on_unsupported_tool:
                    warn(f'{tool} is not a supported ARG annotation tool. Skipping this row.')
                    input_genes.append(None)
                    continue
                else:
                    sys.stderr.write(f'{analysis_software} is not a supported ARG annotation tool\n')
                    sys.stderr.write(f'argNorm can only map genes from the following tools: {list(input_id_lookup.keys())}\n')
                    sys.exit(1)

            if analysis_software == 'groot':
                if '~~~' in row['gene_name']:
                    input_genes.append(input_id_lookup[analysis_software]['groot-argannot'](row))
                elif '.' in row['gene_name']:
                    input_genes.append(input_id_lookup[analysis_software]['groot-card'](row))
                elif 'groot-db_' in row['gene_name']:
                    input_genes.append(input_id_lookup[analysis_software]['groot-db'](row['gene_name']))
                else:
                    input_genes.append(input_id_lookup[analysis_software]['groot-resfinder'](row))
            elif analysis_software == 'abricate':
                try:
                    database = row['reference_database_id']
                except KeyError:
                    try:
                        database = row['reference_database_name']
                    except KeyError:
                        sys.stderr.write(f'An unrecognized hamronization format has been detected. Please use hamronization v1.1.8 or v1.0.4')
                        sys.exit(1)
                input_genes.append(input_id_lookup[analysis_software][database](row))
            else:
                input_genes.append(input_id_lookup[analysis_software](row))
        return pd.Series(input_genes)

    def load_mapping_table(self):
        mapping_tables = []
        for db in DATABASES:
            table = get_aro_mapping_table(db)

            if db == 'resfinder':
                table.index = ResFinderNormalizer.preprocess_ref_genes(table.index)
            elif db == 'ncbi':
                table.index = table.index.str.split('|').str[2]
                protein_id_table = get_aro_mapping_table(db)
                protein_id_table.index = protein_id_table.index.str.split('|').str[1]
                table = pd.concat([table, protein_id_table])
            elif db == 'groot':
                table.index = table.index.str.split('|').str[0]
            elif db == 'megares':
                table.index = table.index.str.split('|').str[0]
            elif db == 'argannot':
                table.index = table.index.map(AbricateNormalizer.preprocess_argannot_ref_gene)
            elif db == 'resfinderfg':
                table.index = table.index.str.split('|').str[1]

            mapping_tables.append(table)

        mapping_table = pd.concat(mapping_tables)
        mapping_table = mapping_table[~mapping_table.index.duplicated(keep='last')]
        return mapping_table


