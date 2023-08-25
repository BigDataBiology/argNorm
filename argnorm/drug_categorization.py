# These functions are run in the run function of the BaseNormalizer.
# After ARO numbers are obtained, these functions can be executed - hence should be independent of which db is used.

import pronto

# Load the ArgNorm ontology from the 'aro.obo' file
ARO = pronto.Ontology.from_obo_library('aro.obo')

def get_immediate_drug_classes(aro_num: str) -> list[tuple]:
    '''
    Description: Gets the drug classes to which a gene confers resistance to.
    Only lists the drug class column in the CARD db.

    Parameters:
        aro_num (str): ARO number. Needs to be in the form 'ARO:number'.

    Returns:
        drug_classes_list (list[tuple]): 
            A two-dimensional list where each inner list represents a drug class.
            Each inner list contains the ARO number and name of the drug class in that order. [ARO:number, name].
    '''

    gene = ARO[aro_num]

    confers_resistance_to_drug_class = any(r.name == 'confers_resistance_to_drug_class' for r in gene.relationships)
    confers_resistance_to_antibiotic = any(r.name == 'confers_resistance_to_antibiotic' for r in gene.relationships)

    drug_classes = []

    if confers_resistance_to_drug_class:
        for drug_class in gene.relationships[ARO.get_relationship('confers_resistance_to_drug_class')]:
            drug_classes.append([drug_class.id, drug_class.name])
    
    if confers_resistance_to_antibiotic:
        for drug_class in gene.relationships[ARO.get_relationship('confers_resistance_to_antibiotic')]:
            drug_classes.append([drug_class.id, drug_class.name])
    
    return drug_classes

def get_drug_class_category(drug_classes_list: list[tuple]) -> list[str]:
    '''
    Description: Gives a list of categories of drug classes, e.g. cephem and penam are categorized as beta_lactam antibiotics.

    Parameters:
        drug_classes_list (list[tuple]): 
            A two-dimensional list where each inner list represents a drug class.
            Each inner list contains the ARO number and name of the drug class in that order. [ARO:number, name].
            Designed to use the return value of the function 'get_immediate_drug_classes'.
    
    Returns:
        drug_class_categories (list[str]):
            A list containing the names of the drug class categories of each drug class given as input to the function.
            Order of the names of the drug class categories corresponds to the order in which the drug classes were given
            to the function in the drug_classes_list.
    '''
    drug_class_categories = []

    for drug_class in drug_classes_list:
        drug_class_instance = ARO[drug_class[0]]
        drug_class_instance_superclasses = list(drug_class_instance.superclasses())
        superclasses_len = len(drug_class_instance_superclasses)

        if superclasses_len >= 3:
            drug_class_categories.append(drug_class_instance_superclasses[superclasses_len - 3].name)
        else:
            drug_class_categories.append(drug_class_instance_superclasses[0].name)
        
    return drug_class_categories
