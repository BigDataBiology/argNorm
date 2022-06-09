import pandas as pd


aros = pd.read_csv('ArgNorm/data/sarg_ARO_mapping.tsv', sep='\t')
ref = pd.read_csv('../quick_amr_db_harmonisation/Ublastx_stageone2.3/DB/structure_20181107.LIST', sep='\t')
ref.Corresponding_ids = ref.Corresponding_ids.apply(eval)
mapping = dict(zip(aros['Original ID'], 'ARO:'  + aros.ARO.astype(str).apply(lambda x: x.split('.')[0])))
ref['ARO'] = ref.Corresponding_ids.apply(lambda x: {mapping[i] for i in x})
ref.to_csv('sarg.orfs_ARO_mapping.tsv', sep='\t')