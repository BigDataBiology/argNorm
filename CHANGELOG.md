# Changelog

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

## Unreleased

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
- A file called `lib.py` will be introduced so that users can use argNorm as a library more easily.
- Users can import the `map_to_aro` function using `from argnorm.lib import map_to_aro`. The function takes a gene name as input, maps the gene to the ARO and returns a pronto term object with the ARO mapping.
- The `get_aro_mapping_table` function, previously within the BaseNormalizer class, has also been moved to `lib.py` to give users the ability to access the mapping tables being used for normalization.
- With the introduction of `lib.py`, users will be able to access core mapping utilities through `argnorm.lib`, drug categorization through `argnorm.drug_categorization`, and the traditional normalizers through `argnorm.normalizers`.
