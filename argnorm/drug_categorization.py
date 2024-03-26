import pronto
from typing import List

# Load the ARO ontology from internet
ARO = pronto.Ontology.from_obo_library('aro.obo')

def confers_resistance_to(aro_num: str) -> List[str]:
    '''
    Description: Returns a list of the drugs/antibiotics to which a gene confers resistance to.

    Parameters:
        aro_num (str): ARO number. Needs to be in the form 'ARO:number'.

    Returns:
        drugs_list (list[str]):
            A list with ARO number of the drugs/antibiotics to which the input gene confers resistance to.
    '''

    if aro_num not in ARO.terms():
        return []
    
    if type(aro_num) == float or type(aro_num) == int:
        aro_num = 'ARO:' + str(aro_num)

    drugs_list = []
    found_drugs_list = False
    gene = ARO[aro_num]

    while not found_drugs_list:
        confers_resistance_to_drug_class = any(r.name == 'confers_resistance_to_drug_class' for r in gene.relationships)
        confers_resistance_to_antibiotic = any(r.name == 'confers_resistance_to_antibiotic' for r in gene.relationships)

        if confers_resistance_to_drug_class:
            found_drugs_list = True
            for drug in gene.relationships[ARO.get_relationship('confers_resistance_to_drug_class')]:
                drugs_list.append(drug.id)

        if confers_resistance_to_antibiotic:
            found_drugs_list = True
            for drug in gene.relationships[ARO.get_relationship('confers_resistance_to_antibiotic')]:
                drugs_list.append(drug.id)

        # ARO:1000001 is 'process or component of antibiotic biology or chemistry'
        if gene.id == 'ARO:1000001':
            break

        gene = list(gene.superclasses(1))[1]
    
    return sorted(drugs_list)

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