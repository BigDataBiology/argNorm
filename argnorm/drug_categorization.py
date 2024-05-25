from typing import List
from . import lib

ARO = lib.get_aro_ontology()
confers_resistance_to_drug_class_rel = ARO.get_relationship('confers_resistance_to_drug_class')
confers_resistance_to_antibiotic_rel = ARO.get_relationship('confers_resistance_to_antibiotic')

def confers_resistance_to(aro_num: str) -> List[str]:
    '''
    Description: Returns a list of the drugs/antibiotics to which a gene confers resistance to.

    Parameters:
        aro_num (str): ARO number. Needs to be in the form 'ARO:number'.

    Returns:
        target (list[str]):
            A list with ARO number of the drugs/antibiotics to which the input gene confers resistance to.
    '''

    target = set()

    for term in ARO[aro_num].superclasses():
        for drug in term.relationships.get(confers_resistance_to_drug_class_rel, []):
            target.add(drug.id)

        for drug in term.relationships.get(confers_resistance_to_antibiotic_rel, []):
            target.add(drug.id)

        if target:
            break

    return sorted(target)

def drugs_to_drug_classes(drugs_list: List[str]) -> List[str]:
    '''
    Description: Returns a list of categories of drug classes, e.g. cephem and penam are categorized as beta_lactam antibiotics.

    Parameters:
        drugs_list (list[str]):
            A list with ARO number of the drugs/antibiotics to which a gene confers resistance to.

    Returns:
        drug_classes (list[str]):
            A list containing the ARO numbers of the drug classes of each drug given as input to the function.
            Order of the ARO numbers of the drug classes corresponds to the order in which the drugs were given
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
