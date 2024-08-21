# Changelog

## Unreleased

## 0.6.0 - 21 August 2024

#### GROOT support
- argNorm supports the GROOT v1.1.2 ARG annotation tool: https://github.com/will-rowe/groot
- GROOT support is via the `GrootNormalizer` (for use in python scripts) and the `groot` tool parameter with the `groot-db`, `groot-core-db`, `groot-argannot`, `groot-card`, and `groot-resfinder` `db` parameters in the CLI.

#### Other
- `__version__` attribute added to the package (accessible as `argnorm.__version__` or `argnorm.lib.__version__`)
- Use atomic writing for outputs (https://github.com/untitaker/python-atomicwrites/tree/master)


### funcscan integration
- argNorm has been included as an nf-core module: https://nf-co.re/modules/argnorm/
- argNorm will also be available on the funcscan pipeline: https://github.com/nf-core/funcscan/pull/410

### DB harmonisation
- SARG db link was changed in `crude_db_harmonisation` to https://raw.githubusercontent.com/xinehc/args_oap/a3e5cff4a6c09f81e4834cfd9a31e6ce7d678d71/src/args_oap/db/sarg.fasta as old link (Galaxy instance, http://smile.hku.hk/SARGs) is down
- RGI outputs in `crude_db_harmonisation` are concatenated so frequencies of `perfect`, `strict`, and `loose` hits can be calculated from concatenated file 

## 0.5.0 - 26 June 2024

### Update drug categorization
- confers_resistance_to() now gets drugs for the whole AMR gene family. For example, OXA-19 previously only returned cephalosporin and penam, but now will also return oxacillin (from AMR gene family).
- Implementation of drugs_to_drug_classes() has also been fixed. Previously, the drug class was obtained from the superclasses of the drugs list passed without a thorough check if the drug class was the immediate child of 'antibiotic molecule'. These checks have now been put in place.
- drugs_to_drug_classes() also uses the 'has_part' ARO relationship now to get drug classes for antibiotic mixtures. In case of antibiotic mixtures, the drug classes of the drugs associated with 'has_part' are returned rather than 'antibiotic mixture' (ARO:3000707).
- 'antibiotic mixture' will not be reported as a drug class, rather the individual antibiotic classes making up the antibiotic mixture will be reported.

### Manual curation
- argannot_curation: (Tet)tetH:EF460464:6286-7839:1554 was incorrectly annotated as ARO:3004797 which is a beta-lactamase due to a loose RGI hit. This was manually curated to ARO:3000175.
- deeparg, megares, resfinderfg & sarg curation: ARO:3004445 -> ARO:3005440, this was due to a change in the ARO and the ARO number for the RSA2 gene changing. **db_harmonisation must change to take this into account**

#### Incorrectly curated genes.
- Previously, these were directly mapped to drug classes. Correct parent ARO term has now been given.
- resfinder_curation: grdA_1_QJX10702 -> 3007380 & EstDL136_1_JN242251 -> 3000557
- megares_curation: MEG_2865|Drugs|Phenicol|Chloramphenicol_hydrolase|ESTD -> 3000557

### Handle AROs as string rather than int in get_aro_mapping_table()

AROs were previously handled as 'int' in the get_aro_mapping_table() function and this posed challenges when ARO numbers such as 'ARO:0010004' (leading zeros are cut). To fix this, AROs are now treated as strings so leading zeros can be maintained.

## 0.4.0 - 10 June 2024

- Bundle a specific version of ARO with the package instead of downloading it from the internet (ensures reproducibility)
- Add missing ARO mappings to manual curation.

    Additions:
    - argannot_curation: (Phe)cpt_strepv:U09991:AAB36569:1412-1948:537 -> ARO: 3000249
    - megares_curation: MEG_2114 -> ARO: 3000249, MEG_2430 -> ARO: 3000016, MEG_985 -> ARO: 3000229, MEG_2865 -> ARO:3000387, MEG_7974 -> ARO:3000076
    - sarg_curation: AM180355.1.gene2260.p01 -> ARO: 3000250
    - resfinder_curation: dldHA2X_1_AL939117 -> ARO:3003970, grdA_1_QJX10702 -> ARO:3007382, EstDL136_1_JN242251 -> ARO:3000387
    - resfinderfg_curation: UDP-N-acetylmuramoyl-tripeptide--D-alanyl-D-alanine ligase|KF629588.1|pediatric_fecal_sample|CYC -> ARO:3003970
- Command line tool accept database/tool names in case-independent way (by @sebastianLedzianowski)
- `lib.map_to_aro` returns `None` if there is no mapping (raises an exception if the name is missing)


## 0.3.0 - 27 April 2024

### Handling gene clusters & reverse complements in resfinder
- Resfinder has gene clusters which can't be passed through RGI using 'contig' mode.
- Gene clusters were identified and were manually assigned ARO numbers.
- A seperate file with manual curation for gene clusters and RCs was created, and their AROs were updated after concatenating RGI results and genes not in RGI results.
- 40 gene clusters present.
- 9 genes in reverse complement form also present.
- RC genes were manually curated.

### Using amino acid file for argannot & resfinder rather than nucleotide file
- ARG-ANNOT and Resfinder are comprised of coding sequences. The data wasn't being handled properly before as contig mode was used when passing coding sequences to RGI. Now, the amino acid versions of ARG-ANNOT & Resfinder are used with protein mode when running the database in RGI.
- ARG-ANNOT AA file is available online. Resfinder AA file is generated using biopython.
- One to many ARO mapping such as NG_047831:101-955 to Erm(K) and almG in ARG-ANNOT eliminated as protein mode used
- A total of 10 ARO mappings changed in ARG-ANNOT

### argnorm.lib: Making argNorm more usable as a library
- Introduce `argnorm.lib` module
- Users can import the `map_to_aro` function from `argnorm.lib`. The function takes a gene name as input, maps the gene to the ARO and returns a pronto term object with the ARO mapping.
- The `get_aro_mapping_table` function, previously within the BaseNormalizer class, has also been moved to `lib.py` to give users the ability to access the mapping tables being used for normalization.
- With the introduction of `lib.py`, users will be able to access core mapping utilities through `argnorm.lib`, drug categorization through `argnorm.drug_categorization`, and the traditional normalizers through `argnorm.normalizers`.


## 0.2.0 - 26 March 2024

#### ARO Mapping & Normalization

- Updated mappings and manual curation tables for latest RGI
- Hamronized ResFinderFG support
- Removed python syntax in output

#### Drug Categorization

- Improved drug categorization by using superclasses whenever direct drug categorization is not possible
- Added better column headings for drug categorization (confers_resistance_to and resistance_to_drug_class)

#### Testing

- Improved pytest testing
- Added integration tests

## 0.1.0 - 20 December, 2023

- Added hamronized support for AMRFinderPlus
- Fixed ARO:nan issue (added manually curated mapping tables and integrated it with normalizers)
- Added drug categorization feature and integrated it with normalizers
- Added AMRFinderPlusNormalizer, ResFinderNormalizer
- Added specific smoke tests for ARGSOAPNormalizer, DeepARGNormalizer, AbricateNormalizer, AMRFinderPlusNormalizer and ResFinderNormalizer

## 0.0.1 - 13 June, 2022

- First release
- Initial source code started
- Normalizers: added BaseNormalizer, ARGSOAPNormalizer, DeepARGNormalizer, AbricateNormalizer
- Testing: added basic ARO column test
