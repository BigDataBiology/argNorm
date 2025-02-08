from typing import List
from . import lib
import pronto

ARO = lib.get_aro_ontology()

confers_resistance_to_drug_class_rel = ARO.get_relationship('confers_resistance_to_drug_class')
confers_resistance_to_antibiotic_rel = ARO.get_relationship('confers_resistance_to_antibiotic')
has_part_rel = ARO.get_relationship('has_part')
part_of_rel = ARO.get_relationship('part_of')
regulates_rel = ARO.get_relationship('regulates')
participates_in_rel = ARO.get_relationship('participates_in')

_antibiotic_class_ids = frozenset(ar.id for ar in ARO['ARO:1000003'].subclasses(distance=1, with_self=False))
def _is_drug_class(ar: str) -> bool:
    return ar in _antibiotic_class_ids


def _get_drug_classes(super_classes_list: List[str]) -> List[pronto.Term]:
    """
    - Helper function to traverse up and record immediate child of 'antibiotic molecule' in ARO
    - Traverses up ARO until immediate child of 'antibiotic molecule' class reached and 'antibiotic mixture' class not reached
    - antibiotic molecule -> ARO:1000003
    - antibiotic mixture -> ARO:3000707
    """
    return [ar for ar in super_classes_list if ar.id in _antibiotic_class_ids and ar.id != 'ARO:3000707']


def _get_drugs(aro_num: str) -> List[pronto.Term]:
    '''
    Description: Returns a list of the drugs/antibiotics to which a gene confers resistance to.

    Parameters:
        aro_num (str): ARO number. Needs to be in the form 'ARO:number'.

    Returns:
        target (list[pronto.Term]):
            ARO nodes the drugs/antibiotics to which the input gene confers resistance to.
    '''

    target = set()

    for superclass in ARO[aro_num].superclasses():
        for drug in superclass.relationships.get(confers_resistance_to_drug_class_rel, []):
            target.add(drug)

        for drug in superclass.relationships.get(confers_resistance_to_antibiotic_rel, []):
            target.add(drug)

        for rel in [regulates_rel, participates_in_rel, part_of_rel]:
            for term in superclass.relationships.get(rel, []):
                target.update(_get_drugs(term.id))

    return sorted(target)


def confers_resistance_to(aro_num: str) -> List[str]:
    # some gene superclasses can map to drugs which are immediate children of 'antibiotic molecule'
    # only use these if they are not redundant with other drugs

    drugs = set(_get_drugs(aro_num))
    drug_classes = set(d for d in drugs if _is_drug_class(d.id))
    drugs = drugs - drug_classes

    for drug in drugs:
        cur_classes = _get_drug_classes(drug.superclasses(with_self=False))
        drug_classes -= set(cur_classes)

    drugs.update(drug_classes)

    return sorted(d.id for d in drugs)


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
        drug = ARO[drug]
        cur_drug_classes = _get_drug_classes(drug.superclasses())

        for part in drug.relationships.get(has_part_rel, []):
            cur_drug_classes += _get_drug_classes(part.superclasses())
        drug_classes.extend(set(cur_drug_classes))

    return sorted(d.id for d in drug_classes)

