# These functions are run in the run function of the BaseNormalizer.
# After ARO numbers are obtained, these functions can be executed - hence should be independent of which db is used.

import pronto
from typing import List, Tuple

# Load the ARO ontology from internet
ARO = pronto.Ontology.from_obo_library('aro.obo')

def confers_resistance_to(aro_num: str) -> List[Tuple]:
    '''
    Description: Gets the drugs/antibiotics to which a gene confers resistance to.
    Only lists the drug class column in the CARD db.

    Parameters:
        aro_num (str): ARO number. Needs to be in the form 'ARO:number'.

    Returns:
        drugs_list (list[tuple]):
            A two-dimensional list where each inner list represents a drug.
            Each inner list contains the ARO number and name of the drug in that order. [ARO:number, name].
    '''

    # Some databases don't provide aro numbers as strings.
    # Converting those aro numbers to pronto's desired format.
    if type(aro_num) == float or type(aro_num) == int:
        aro_num = 'ARO:' + str(aro_num)

    # If dealing with aro nans, the final drug class categorization will give [].
    # Hence, immediate drug classes categorization also gives empty list.
    if aro_num == 'ARO:nan':
        return []

    # Returning empty list if aro number not in ARO ontology.
    if aro_num not in ARO.terms():
        return []

    gene = ARO[aro_num]

    confers_resistance_to_drug_class = any(r.name == 'confers_resistance_to_drug_class' for r in gene.relationships)
    confers_resistance_to_antibiotic = any(r.name == 'confers_resistance_to_antibiotic' for r in gene.relationships)

    drugs_list = []

    if confers_resistance_to_drug_class:
        for drug in gene.relationships[ARO.get_relationship('confers_resistance_to_drug_class')]:
            drugs_list.append(drug.id)

    if confers_resistance_to_antibiotic:
        for drug in gene.relationships[ARO.get_relationship('confers_resistance_to_antibiotic')]:
            drugs_list.append(drug.id)

    return sorted(drugs_list)

def drugs_to_drug_classes(drugs_list: List[Tuple]) -> List[Tuple]:
    '''
    Description: Gives a list of categories of drug classes, e.g. cephem and penam are categorized as beta_lactam antibiotics.

    Parameters:
        drugs_list (list[tuple]):
            A two-dimensional list where each inner list represents a drug.
            Each inner list contains the ARO number and name of the drug in that order. [ARO:number, name].
            Designed to use the return value of the function 'confers_resistance_to'.

    Returns:
        drug_classes (list[str]):
            A list containing the names of the drug classes of each drug given as input to the function.
            Order of the names of the drug classes corresponds to the order in which the drugs were given
            to the function in the drugs_list.
    '''
    drug_classes = []

    for drug in drugs_list:
        drug_instance = ARO[drug]
        drug_instance_superclasses = list(drug_instance.superclasses())
        superclasses_len = len(drug_instance_superclasses)

        if superclasses_len >= 3:
            drug_classes.append(drug_instance_superclasses[superclasses_len - 3].id)
        else:
            drug_classes.append(drug_instance_superclasses[0].id)

    return sorted(drug_classes)
