from .normalizers import ARGSOAPNormalizer, \
    DeepARGNormalizer, AbricateNormalizer, ResFinderNormalizer, AMRFinderPlusNormalizer

def normalize(ifile, tool : str, database : str, is_hamronized : bool, mode : str):
    '''Normalize ARG tables

    Parameters
    ----------

    ifile : input file
    tool : str
    database : str
        Database to use
    is_hamronized : bool
        Whether input has been run through hAMRonization already
    mode : str
    '''
    normalizer = {
        'abricate': AbricateNormalizer,
        'amrfinderplus': AMRFinderPlusNormalizer,
        'argsoap': ARGSOAPNormalizer,
        'deeparg': DeepARGNormalizer,
        'resfinder': ResFinderNormalizer,
    }.get(tool)

    if normalizer is None:
        raise ValueError('Please specify a correct tool name.')
    norm = normalizer(database=database, is_hamronized=is_hamronized, mode=mode)
    return norm.run(ifile)

