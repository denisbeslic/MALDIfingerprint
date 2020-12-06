# Identification of Antibody Subclasses

To determine the subclass of an antibody, the mass spectrum is compared with a library of simulated mass spectra of each antibody class. 

### Input
- directory with experimental MS files in .txt format with mass and intensity columns
- directory with theoretical MS files in .txt format with mass and intensity columns

### Output
- overview table showing number of peaks that between all experimental an theoretical MS files in .txt format
- experimental MS files with additional "subclass" column 

## Generate library files
1. Find sequence of subclass on Uniprot or IMGT
2. Use MS-Digest to create theoretical spectra
3. 