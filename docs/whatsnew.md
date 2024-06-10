# What's New

## 0.4.0 - 10 June

### User-facing changes
- Command line tool accept database/tool names in case-independent way (by @sebastianLedzianowski)
- A few ARO mappings were missing in the manual curation and they have been added. See the [CHANGELOG](https://github.com/BigDataBiology/argNorm/blob/main/CHANGELOG.md) for more details.
- `lib.map_to_aro` returns `None` if there is no mapping (raises an exception if the name is missing)

### Internal changes

- Bundle a specific version of ARO with the package instead of downloading it from the internet (ensures reproducibility)

## 0.3.0 - 27 April 2024

### User-facing changes
- Improved Resfinder mappings (40 gene clusters and 9 reverse complements were manually curated)
- Updated ARG-ANNOT mappings (a total of 10 mappings changed)
- argNorm is now more usable as a library
- Remove warnings

### Internal changes
- Code has been refactored to be simpler


## 0.2.0 - 26 March 2024

#### ARO Mapping & Normalization

- Updated mappings and manual curation tables for latest RGI
- Hamronized ResFinderFG support
- Removed python syntax in output

#### Drug Categorization

- Improved drug categorization by using superclasses whenever direct drug categorization is not possible
- Added better column headings for drug categorization (confers_resistance_to and resistance_to_drug_class)

#### Internal Improvements: Testing

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
