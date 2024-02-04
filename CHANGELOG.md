# Changelog

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
- Changed manual curation tables to use parent term ARO if a direct mapping can't be done
- Hamronized ResFinderFG support
- New example outputs created using latest argNorm version from examples directory
- Added ARO numbers to OVERALL CATEGORY OF DRUG CLASS