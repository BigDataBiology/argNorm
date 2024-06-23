from typing import List
from . import lib

ARO = lib.get_aro_ontology()
confers_resistance_to_drug_class_rel = ARO.get_relationship('confers_resistance_to_drug_class')
confers_resistance_to_antibiotic_rel = ARO.get_relationship('confers_resistance_to_antibiotic')
has_part_rel = ARO.get_relationship('has_part')

def _get_drug_classes(super_classes_list: List[str]) -> List[str]:
    """
    - Helper function to traverse up and record immediate child of 'antibiotic molecule' in ARO
    - Traverses up ARO until immediate child of 'antibiotic molecule' class reached and 'antibiotic mixture' class not reached
    - antibiotic molecule -> ARO:1000003
    - antibiotic mixture -> ARO:3000707
    """
    output = []

    for super_class in super_classes_list:
        super_class_classes = list(super_class.superclasses(1))
        antibiotic_molecule_node = [ARO['ARO:1000003']]

        # checking if immediate child of 'antibiotic molecule' is reached & it is not 'antibiotic mixture'
        if super_class_classes[1:] == antibiotic_molecule_node and super_class.id != 'ARO:3000707':
            output.append(super_class.id)

    return output

def confers_resistance_to(aro_num: str) -> List[str]:
    '''
    Description: Returns a list of the drugs/antibiotics to which a gene confers resistance to.

    Parameters:
        aro_num (str): ARO number. Needs to be in the form 'ARO:number'.

    Returns:
        target (list[str]):
            A list with ARO number of the drugs/antibiotics to which the input gene confers resistance to.
    '''
    antibiotic_molecule_node = [ARO['ARO:1000003'], ARO['ARO:1000001']]
    # some gene superclasses can map to drugs which are immediate children of 'antibiotic molecule'
    # only use these if no other drugs can be found, as this information will be present in
    # drugs to drug classes
    backup_drugs = []
    target = set()
    for term in ARO[aro_num].superclasses():
        for drug in term.relationships.get(confers_resistance_to_drug_class_rel, []):
            if list(ARO[drug.id].superclasses())[1:] == antibiotic_molecule_node:
                backup_drugs.append(drug.id)
            else:
                target.add(drug.id)

        for drug in term.relationships.get(confers_resistance_to_antibiotic_rel, []):
            if list(ARO[drug.id].superclasses())[1:] == antibiotic_molecule_node:
                backup_drugs.append(drug.id)
            else:
                target.add(drug.id)

    if not target:
        target.update(backup_drugs)

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
    temp_drug_classes = []

    for drug in drugs_list:
        drug_instance = ARO[drug]
        drug_instance_superclasses = list(drug_instance.superclasses())
        temp_drug_classes += _get_drug_classes(drug_instance_superclasses)

        has_part_nodes = drug_instance.relationships.get(has_part_rel, [])
        for has_part_node in has_part_nodes:
            has_part_node_superclasses = list(has_part_node.superclasses())[1:]

            for super_class in has_part_node_superclasses:
                super_class_categories = list(super_class.superclasses())
                temp_drug_classes += _get_drug_classes(super_class_categories)

        if temp_drug_classes == []:
            temp_drug_classes.append(drug_instance.id)

        drug_classes += list(set(temp_drug_classes))
        temp_drug_classes = []

    return sorted(drug_classes)
