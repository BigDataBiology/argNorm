"""
A TEMPORARY bypass from the main normalization pipeline if drug categorization
needs to be implemented on already ARO normalized output.
"""

from drug_categorization import get_immediate_drug_classes
from drug_categorization import get_drug_class_category
import pandas as pd

IMMEDIATE_DRUG_CLASS_COL_HEADING = 'CONFERS RESISTANCE TO IMMEDIATE DRUG CLASS'
DRUG_CLASS_CATEGORY_COL_HEADING = 'OVERALL CATEGORY OF DRUG CLASS'
TARGET_ARO_COL = 'ARO'


def load_input(input_file):
    """
    Customize this when it fails to parse the input data.
    """
    return pd.read_csv(input_file, sep='\t')


def initial_drug_categorization(aro_list):
    immediate_drug_classes = []
    for aro in aro_list:
        immediate_drug_classes.append(get_immediate_drug_classes(aro))

    return immediate_drug_classes


def final_drug_categorization(immediate_drug_classes_col):
    drug_class_catogeries_col = []
    for drug_classes in immediate_drug_classes_col:
        drug_class_catogeries_col.append(get_drug_class_category(drug_classes))

    return drug_class_catogeries_col


def perform_drug_categorization_on_aro_normalized_output(input_file: str):
    normalized_output = load_input(input_file)
    normalized_output[IMMEDIATE_DRUG_CLASS_COL_HEADING] = initial_drug_categorization(normalized_output[TARGET_ARO_COL])
    normalized_output[DRUG_CLASS_CATEGORY_COL_HEADING] = final_drug_categorization(normalized_output[IMMEDIATE_DRUG_CLASS_COL_HEADING])
    return normalized_output
