from .normalizers import ARGSOAPNormalizer, \
    DeepARGNormalizer, AbricateNormalizer, ResFinderNormalizer, AMRFinderPlusNormalizer, GrootNormalizer, HamronizationNormalizer

def normalize(ifile, tool : str, database : str, skip_on_unsupported_tool : bool):
    '''Normalize ARG tables

    Parameters
    ----------

    ifile : input file
    tool : str
    database : str
        Database to use
    is_hamronized : bool
        Whether input has been run through hAMRonization already
    '''
    normalizer = {
        'abricate': AbricateNormalizer,
        'amrfinderplus': AMRFinderPlusNormalizer,
        'argsoap': ARGSOAPNormalizer,
        'deeparg': DeepARGNormalizer,
        'resfinder': ResFinderNormalizer,
        'groot': GrootNormalizer,
        'hamronization': HamronizationNormalizer
    }.get(tool)

    if normalizer is None:
        raise ValueError('Please specify a correct tool name.')
    
    if normalizer == HamronizationNormalizer:
        norm = normalizer(database=database, skip_on_unsupported_tool=skip_on_unsupported_tool)
    else:
        norm = normalizer(database=database)
    
    return norm.run(ifile)

