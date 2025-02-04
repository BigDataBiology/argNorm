import pandas as pd
try:
    pd.options.mode.copy_on_write = True
except pd.errors.OptionError:
    pass
from .drug_categorization import confers_resistance_to, drugs_to_drug_classes
from .lib import get_aro_mapping_table
from .lib import MAPPING_TABLE_ARO_COL, TARGET_ARO_COL, DATABASES
from warnings import warn

# Column headings for drug categorization output
CONFERS_RESISTANCE_TO_COL = 'confers_resistance_to'
RESISTANCE_TO_DRUG_CLASSES_COL = 'resistance_to_drug_classes'

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

        aro_table = get_aro_mapping_table(self.database)
        aro_table.set_index(self.preprocess_ref_genes(
            aro_table.index
        ), inplace=True)
        mapping = aro_table[MAPPING_TABLE_ARO_COL].to_dict()
        original_annot[TARGET_ARO_COL] = input_genes.map(mapping)

        # Drug categorization
        original_annot[CONFERS_RESISTANCE_TO_COL] = original_annot[TARGET_ARO_COL].map(
                lambda a: ','.join(confers_resistance_to(a)),
                na_action='ignore')
        original_annot[RESISTANCE_TO_DRUG_CLASSES_COL] = original_annot[TARGET_ARO_COL].map(
                lambda a: ','.join(drugs_to_drug_classes(confers_resistance_to(a))),
                na_action='ignore')

        return original_annot

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
        database = 'sarg'
        super().__init__(database)

    def get_input_ids(self, itable):
        return itable[1]

    def load_input(self, input_file):
        return pd.read_csv(input_file, sep='\t', header=None)

class DeepARGNormalizer(BaseNormalizer):
    def __init__(self, database=None) -> None:
        database = 'deeparg'
        super().__init__(database)

    def get_input_ids(self, itable):
        return itable['best-hit']

class ResFinderNormalizer(BaseNormalizer):
    def __init__(self, database=None) -> None:
        database = 'resfinder'
        super().__init__(database)

    def get_input_ids(self, itable):
        gene_identifier = 'Resistance gene'
        accession = 'Accession no.'
        return pd.Series(itable[gene_identifier] + '_' + itable[accession])

    def preprocess_ref_genes(self, ref_genes):
        split_genes = ref_genes.str.split('_')
        return pd.Series(split_genes.str[0] + '_' + split_genes.str[-1])

class AMRFinderPlusNormalizer(BaseNormalizer):
    def __init__(self, database=None) -> None:
        database = 'ncbi'
        super().__init__(database)

    def get_input_ids(self, itable):
        gene_identifier = 'Sequence name'
        accession = 'Accession of closest sequence'
        return pd.Series(itable[accession] + '|' + itable[gene_identifier].str.replace(' ', '_'))

    def preprocess_ref_genes(self, ref_genes):
        split_genes = ref_genes.str.split('|')
        return pd.Series(split_genes.str[1] + '|' + split_genes.str[-1])

class AbricateNormalizer(BaseNormalizer):
    def __init__(self, database=None) -> None:
        if database not in ['ncbi', 'deeparg', 'resfinder', 'megares', 'argannot']:
            raise Exception(f'{database} is not a supported database.')

        super().__init__(database)

    def get_input_ids(self, itable):
        col = dict(
            ncbi='PRODUCT',
            deeparg='best-hit',
            resfinder=dict(gene_identifier='GENE', accession='ACCESSION'),
            megares='ACCESSION',
            argannot='ACCESSION'
        )[self.database]

        if self.database == 'resfinderfg':
            return itable[col].str.split('|').str[1]
        
        if self.database == 'resfinder':          
            return pd.Series(itable[col['gene_identifier']] + '_' + itable[col['accession']])
            
        if self.database == 'ncbi':
            return pd.Series(itable[col].str.replace(' ', '_'))
        
        return itable[col]

    def preprocess_argannot_ref_genes(self, ref_gene):
        split_str = ref_gene.split(':')
        if not str(split_str[2][0]).isnumeric() and not '-' in split_str[2]:
            return ':'.join([split_str[1], split_str[3]])
        return ':'.join(ref_gene.split(':')[1:3])

    def preprocess_ref_genes(self, ref_genes):
        process_funcs_by_db = dict(
            ncbi=lambda x: x.split('|')[-1],
            deeparg=lambda x: x,
            resfinder=lambda x: '_'.join([x.split('_')[0], x.split('_')[-1]]),
            sarg=lambda x: x,
            megares=lambda x: x.split('|')[0],
            argannot=self.preprocess_argannot_ref_genes,
            resfinderfg=lambda x: x.split('|')[1]
        )
        return ref_genes.map(process_funcs_by_db[self.database])

class GrootNormalizer(BaseNormalizer):
    def __init__(self, database=None) -> None:
        if database not in ['groot-argannot', 'groot-resfinder', 'groot-db', 'groot-core-db', 'groot-card']:
            raise Exception(f'{database} is not a supported database for groot.')

        super().__init__(database)
            
    def load_input(self, input_file):
        return pd.read_csv(input_file, sep='\t', header=None)
    
    def preprocess_groot_db_inputs(self, gene_name):
        processed_gene_name = str(gene_name).split('__')[1]
        if 'card' in  str(gene_name).split('__')[0].lower():
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
    def __init__(self, database=None):
        super().__init__(database)

        self.input_ids = {
            'argsoap': lambda x: x['reference_accession'],
            'deeparg': lambda x: x['gene_name'],
            'resfinder': lambda x: x['gene_name'] + '_' + x['reference_accession'],
            'amrfinderplus': lambda x: str(x['gene_name']).replace(' ', '_'),
            'groot': {
                'groot-card': lambda x: x['gene_name'].split('.')[0], 
                'groot-argannot': lambda x: x['gene_name'].split('~~~')[-1],
                'groot-db': self.preprocess_groot_db_inputs,
                'groot-resfinder': lambda x: x['gene_name']
            },
            'abricate': {
                'sarg': lambda x: x['gene_symbol'],
                'megares': lambda x: x['reference_accession'],
                'argannot': lambda x: x['reference_accession'],
                'resfinderfg': lambda x: x['gene_name'].split('|')[1],
                'deeparg': lambda x: x['gene_name'],
                'resfinder': lambda x: x['gene_name'] + '_' + x['reference_accession'],
                'ncbi': lambda x: str(x['gene_name']).replace(' ', '_')
            }
        }

        self.preprocess_ref_genes_funcs = {
            'resfinder': lambda x: x.split('_')[0] + '_' + x.split('_')[-1],
            'ncbi': lambda x: x.split('|')[-1],
            'groot': lambda x: ':'.join(str(x).split(':')[1:3]) if ':' in x else x,
            'deeparg': lambda x: x,
            'sarg': lambda x: x,
            'megares': lambda x: x.split('|')[0],
            'argannot': self.preprocess_argannot_ref_genes,
            'resfinderfg': lambda x: x.split('|')[1],
        }

    def preprocess_groot_db_inputs(self, gene_name):
        processed_gene_name = str(gene_name).split('__')[1]
        if 'card' in  str(gene_name).split('__')[0].lower():
            processed_gene_name = processed_gene_name.split('|')[-1]
        return processed_gene_name
    
    def preprocess_argannot_ref_genes(self, ref_gene):
        split_str = ref_gene.split(':')
        if not str(split_str[2][0]).isnumeric() and not '-' in split_str[2]:
            return ':'.join([split_str[1], split_str[3]])
        return ':'.join(ref_gene.split(':')[1:3])

    def run(self, input_file : str):
        original_annot = self.load_input(input_file)

        input_genes = []
        for i in original_annot.index:
            analysis_software = str(original_annot.iloc[i]['analysis_software_name']).lower().replace('-', '').replace('_', '').replace(' ', '')
            
            for tool in list(self.input_ids.keys()):
                if tool in analysis_software:
                    analysis_software = tool
                    break

                if 'ncbi' in analysis_software:
                    analysis_software = 'amrfinderplus'
                    break
            else:
                warn(f'{analysis_software} is not a supported ARG annotation tool')
                input_genes.append(None)
                continue

            if analysis_software == 'groot':
                if '~~~' in original_annot.iloc[i]['gene_name']:
                    input_genes.append(self.input_ids[analysis_software]['groot-argannot'](original_annot.iloc[i]))
                elif '.' in original_annot.iloc[i]['gene_name']:
                    input_genes.append(self.input_ids[analysis_software]['groot-card'](original_annot.iloc[i]))
                elif 'groot-db_' in original_annot.iloc[i]['gene_name']:
                    input_genes.append(self.input_ids[analysis_software]['groot-db'](original_annot.iloc[i]['gene_name']))
                else:
                    input_genes.append(self.input_ids[analysis_software]['groot-resfinder'](original_annot.iloc[i]))
            elif analysis_software == 'abricate':
                database = original_annot.iloc[i]['reference_database_id']
                input_genes.append(self.input_ids[analysis_software][database](original_annot.iloc[i]))
            else:
                input_genes.append(self.input_ids[analysis_software](original_annot.iloc[i]))

        mapping_tables = []
        for db in DATABASES:
            table = get_aro_mapping_table(db)
            if db in list(self.preprocess_ref_genes_funcs.keys()):
                table.index = list(map(self.preprocess_ref_genes_funcs[db], table.index))
            mapping_tables.append(table)

        aro_table = pd.concat(mapping_tables)
        mapping = aro_table[MAPPING_TABLE_ARO_COL].to_dict()
        original_annot[TARGET_ARO_COL] = pd.Series(input_genes).map(mapping)

        original_annot[CONFERS_RESISTANCE_TO_COL] = original_annot[TARGET_ARO_COL].map(
                lambda a: ','.join(confers_resistance_to(a)), na_action='ignore')
        original_annot[RESISTANCE_TO_DRUG_CLASSES_COL] = original_annot[TARGET_ARO_COL].map(
                lambda a: ','.join(drugs_to_drug_classes(confers_resistance_to(a))), na_action='ignore')

        return original_annot
