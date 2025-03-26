# Changelog

## Unreleased version

### User-facing changes

- Added columns to output with English names for AROs, drugs, and drug classes (`ARO_name`, `confers_resistance_to_names`, `resistance_to_drug_classes_names`)
- Added `--version` argument

## Version 0.7.0 - 5 March 2025

### User-facing changes

#### Support amrfinderplus v4 alongisde amrfinderplus v3
- Both the latest version of amrfinderplus (v4.0.19) and the previous v3.10.30 are now supported.

#### Added `Cut_Off` column to outputs
- A new column called `Cut_Off` is now present in all argNorm outputs.
- These cut off scores are procured from RGI outputs after running input databases through RGI, and these cut-off scores correspond to RGI's Discovery paradigm (learn more here: https://github.com/arpcard/rgi/blob/master/docs/rgi_main.rst).
- RGI's cut-off scores are of three types: `Loose`, `Strict`, and `Perfect`.
- argNorm also adds another score, `Manual`, to signify genes that are mapped to ARO terms manually without using RGI.


#### Added `HamronizationNormalizer` and `--hamronization_skip_unsupported_tool`
#### Added `HamronizationNormalizer` and `--hamronization_skip_unsupported_tool`
- All hamronized results now go through the `HamronizationNormalizer` class.
- Removed the `is_hamronized` property for all normalizers and removed `--hamronized` flag for CLI.
- HamronizationNormalizer reads a hamronized file line by line, procures input genes, and loads all ARO mapping tables to support hamronized results that combine the outputs from multiple tools and databases.
- For CLI hamronization commands will look like:
```bash
argnorm hamronization -i PATH_TO_INPUT -o PATH_TO_OUTPUT [--hamronization_skip_unsupported_tool]
```
- Combined hamronization results can have ARGs detected by unsupported tools (e.g. staramr). By default, argNorm throws an exception as these are unsupported tools, however, `--hamronization_skip_unsupported_tool` allows users to skip rows with unsupported tools. A warning will be raised rather than an exception.

#### Update `confers_resistance_to()` to use `regulates`, `part_of`, and `participates_in` ARO relationships
Previously, argNorm used the `is_a` ARO relationship along with `confers_resistance_to_drug_class` and `confers_resistance_to_antibiotic` to map ARGs to the drugs they confer resistance to. While this worked well for most genes, some ARGs such as those coding for efflux pumps/proteins (e.g. `ARO:3003548`, `ARO:3000826`, `ARO:3003066`) were previously not mapped to any drugs. This is because none of their superclasses mapped to drugs/antibiotics via `confers_resistance_to_antibiotic` or `confers_resistance_to_drug_class`. However, these genes were related to other ARGs that did map to drugs via the `regulates`, `part_of`, or `participates_in` ARO relationships. argNorm now also utilizes these three relationships to ensure that even if the superclasses (derived using `is_a`) of an ARG don't map to a drug, the gene can be assigned a drug mapping.

#### argNorm mappings updated for CARD and ARO v4.0.0
On 18/12/2024, CARD and ARO were updated to v4.0.0. This brought significant changes to both CARD and the ARO, with 1200+ new beta-lactamase genes being added. argNorm's ARO mappings have now been updated to support CARD & ARO v4.0.0.

> **NOTE: ARO mappings for many ARGs have changed. An ARO mapping from a previous version of argNorm might not be present in the latest ARO or can even point to a completely different gene.**

#### argNorm version in outputs from CLI
The version of argNorm used for normalization is now added to the top of argNorm tsv outputs when generated using the CLI.  For example,

~~~
# argNorm version: 0.6.0
input_sequence_id	input_file_name	gene_symbol	gene_name	reference_database_id ...
. REST OF ARG ANNOTATION OUTPUT TSV TABLE
.
.
~~~

> **NOTE: This is only if argNorm is used on the CLI, if you used argNorm's normalizers there will be no comment with the argNorm version**

> **NOTE: THIS WILL BREAK ANY PREVIOUS SCRIPTS ANALYZING ARGNORM CLI OUTPUTS BEFORE THIS UPDATE! PLEASE USE THE `skiprows=1` ARGUMENT WHEN LOADING ARGNORM OUTPUT DATAFRAMES TO IGNORE THE COMMENT WITH THE ARGNORM VERSION AS SHOWN BELOW:**
> ```
> import pandas as pd
> df = pd.read_csv(<PATH TO ARGNORM OUTPUT>, sep='\t', skiprows=1)
> ```

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
