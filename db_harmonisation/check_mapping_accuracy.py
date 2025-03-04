from argnorm.drug_categorization import confers_resistance_to, drugs_to_drug_classes, ARO
import pandas as pd


def preprocess_mappings_for_tests(mapping, database):
    rgi_mapping = mapping

    # Add drug class information from argnorm
    drug_classes_list = []

    for i in range(rgi_mapping.shape[0]):
        aro = "ARO:" + str(rgi_mapping.iloc[i]['ARO'])
        drugs = confers_resistance_to(aro)
        drug_classes = drugs_to_drug_classes(drugs)
        drug_classes_list.append(list(map(lambda x: ARO[x].name, drug_classes)))


    rgi_mapping['Drug Classes'] = drug_classes_list

    # Add drug class information from original database
    resfinder_antibiotic_classes = pd.read_csv('./drug_classes/resfinder_antibiotic_classes.tsv', sep='\t')
    sarg_antibiotic_classes = pd.read_csv('./drug_classes/SARG_structure.tsv', sep='\t')
    drug_classes = []
    
    for i in range(rgi_mapping.shape[0]):
        if database == 'argannot':
            drug_class = rgi_mapping.iloc[i]['Original ID'].split(')')[0][1:]
        if database == 'ncbi':
            if rgi_mapping.iloc[i]['Original ID'].split('|')[-2]:
                drug_class = rgi_mapping.iloc[i]['Original ID'].split('|')[-2]
            else:
                drug_class = rgi_mapping.iloc[i]['Original ID']
        if database == 'deeparg':
            drug_class = rgi_mapping.iloc[i]['Original ID'].split('|')[-2]
        if database == 'resfinder':
            gene_name = rgi_mapping.iloc[i]['Original ID']
            drug_class = str(resfinder_antibiotic_classes[resfinder_antibiotic_classes['Gene_accession no.'] == gene_name]['Class'].values).replace("['", '').replace("']", '')
        if database == 'resfinder_fg':
            drug_class = rgi_mapping.iloc[i]['Original ID'].split('|')[0]
        if database == 'megares':
            drug_class = rgi_mapping.iloc[i]['Original ID'].split('|')[2]
        if database == 'sarg':
            gene_name = rgi_mapping.iloc[i]['Original ID'].split(' ')[0]
            drug_class = str(sarg_antibiotic_classes[sarg_antibiotic_classes['SARG.Seq.ID'] == gene_name]['Type'].values).replace("['", '').replace("']", '')

        drug_classes.append(drug_class)
    
    rgi_mapping['Original Drug Classes'] = drug_classes
    return rgi_mapping

def check_mapping_accuracy(processed_mappings : pd.DataFrame):
    alternative_ids = {
        "beta-lactam antibiotic": ['Methicillin resistance mecR1 protein', 'Bla', 'beta_lactam', 'beta-lactam', 'Beta-lactamase', 'betalactams', 'beta-lactamase', 'Metallo-beta-lactamase', 'putative peptidoglycan D%2CD-transpeptidase PenA'],
        "aminoglycoside antibiotic": ['aminoglycoside', 'AGly', 'Aminoglycosides', "Streptomycin 3''-adenylyltransferase", 'Gentamicin 3-N-acetyltransferase', 'Bifunctional AAC/APH'],
        "macrolide antibiotic": ['MLS', 'MACROLIDE', 'Macrolide'],
        "tetracycline antibiotic": ['tetracycline', 'Tet', 'TETRACYCLINE', 'Tetracyclines', 'Tetracycline', 'Tetracycline resistance protein', 'Tetracycline repressor protein'],
        "phenicol antibiotic": ['Phe', 'chloramphenicol', 'amphenicol', 'Chloramphenicol acetyltransferase', 'Chloramphenicol acetyltransferase 2', 'florfenicol'],
        "phosphonic acid antibiotic": ['fosfomycin', 'Fosfomycin', 'Fcyn'],
        "sulfonamide antibiotic": ['Sulfonamides', 'Dihydropteroate synthase', 'Folate pathway antagonist'],
        "glycopeptide antibiotic": ['Glycopeptides', 'vancomycin', 'bleomycin', 'D-alanine--D-alanine ligase', 'D-alanine--D-alanine ligase B', 'D-alanine--D-alanine ligase A'],
        "diaminopyrimidine antibiotic": ['Dihydrofolate reductase', 'Trimethoprim', 'Folate pathway antagonist', 'Tmt'],
        "peptide antibiotic": ['bacitracin', 'polymyxin', 'other_peptide_antibiotics', 'COLISTIN', 'lipopeptides', 'COL', 'TUBERACTINOMYCIN', 'edeine', 'defensin', 'Cationic_antimicrobial_peptides'],
        "aminocoumarin antibiotic": ['novobiocin'],
        "nucleoside antibiotic": ['puromycin', 'Nucleosides', 'tunicamycin', 'streptothricin', 'aminoglycoside', 'AGly'],
        "rifamycin antibiotic": ['rifampin'],
        "lincosamide antibiotic": ['lincosamide', 'MLS'],
        "streptogramin antibiotic": ['streptogramin', 'streptogramin A', 'streptogramin B', 'MLS'],
        'nitroimidazole antibiotic': ['Metronidazole', 'Ntmdz'],
        'oxazolidinone antibiotic': ['Oxzln'],
        'fusidane antibiotic': ['fusidic_acid', 'fusaric-acid', 'FUSIDIC_ACID', 'fusidic-acid', 'Fcd'],
        'fluoroquinolone antibiotic': ['Flq', 'Fluoroquinolones', 'PHENICOL/QUINOLONE'],
        'pleuromutilin antibiotic': ['pleuromutilin_tiamulin', 'LINCOSAMIDE/PLEUROMUTILIN'],
    }
    metals = ['mercury_resistance', 'multi-metal_resistance', 'tellurium_resistance', 'tellurium', 'arsenic', 'cadmium', 'copper', 'mercury', 'nickel', 'copper/silver', 'silver', 'cadmium/cobalt/nickel', 'chromate', 'COPPER/GOLD', 'GOLD']
    virulence_genes_or_toxins = ['stx2', 'intimin', 'stx1']
    
    detected_metal_virulence_biocide_genes = []
    drug_class_mismatch_list = []
    
    for i in range(processed_mappings.shape[0]):
        if str(processed_mappings.iloc[i]['Original Drug Classes']).lower() == str(processed_mappings.iloc[i]['Drug Classes']).lower():
            pass
        else:
            if str(processed_mappings.iloc[i]['Original Drug Classes']).lower() in str(processed_mappings.iloc[i]['Drug Classes']).lower():
                pass
            elif 'multidrug' in str(processed_mappings.iloc[i]['Original Drug Classes']).lower() or 'multi-drug_resistance' in str(processed_mappings.iloc[i]['Original Drug Classes']).lower():
                pass
            elif 'macrolide-lincosamide-streptogramin' in str(processed_mappings.iloc[i]['Original Drug Classes']).lower():
                matched = False
                for ii in 'macrolide-lincosamide-streptogramin'.split('-'):
                    if ii in str(processed_mappings.iloc[i]['Drug Classes']).lower():
                        matched = True
                        break
                
                if not matched:
                    drug_class_mismatch_list.append(processed_mappings.iloc[i]['Original ID'])
            elif 'sdia' in str(processed_mappings.iloc[i]['Original ID']).lower() or 'cpxr' in str(processed_mappings.iloc[i]['Original ID']).lower() or 'rosA' in str(processed_mappings.iloc[i]['Original ID']) or 'rosB' in str(processed_mappings.iloc[i]['Original ID']):
                pass
            elif 'penicillin-binding_protein_' in str(processed_mappings.iloc[i]['Original ID']).lower() and 'beta-lactam antibiotic' in str(processed_mappings.iloc[i]['Drug Classes']):
                pass
            elif str(processed_mappings.iloc[i]['Original Drug Classes']).lower() in metals:
                detected_metal_virulence_biocide_genes.append(processed_mappings.iloc[i]['Original ID'])
            elif str(processed_mappings.iloc[i]['Original Drug Classes']).lower() in virulence_genes_or_toxins:
                detected_metal_virulence_biocide_genes.append(processed_mappings.iloc[i]['Original ID'])
            elif str(processed_mappings.iloc[i]['Original Drug Classes']) == "Drug_and_biocide_resistance":
                detected_metal_virulence_biocide_genes.append(processed_mappings.iloc[i]['Original ID'])
            else:
                matched = False
                for id in alternative_ids:
                    if id in processed_mappings.iloc[i]['Drug Classes']:
                        for alternative in alternative_ids[id]:
                            if str(alternative).lower() in str(processed_mappings.iloc[i]['Original Drug Classes']).lower():
                                matched = True
                                break
                        
                if not matched:
                    drug_class_mismatch_list.append(processed_mappings.iloc[i]['Original ID'])
    
    mismatched_genes = processed_mappings.loc[processed_mappings['Original ID'].isin(drug_class_mismatch_list)]
    return {"mismatched_genes": mismatched_genes, "metal_biocide_virulence_genes": detected_metal_virulence_biocide_genes}
