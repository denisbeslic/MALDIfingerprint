# Identification of Antibody Subclasses

To determine the subclass of an antibody, the mass spectrum is compared with a library of simulated mass spectra of each antibody class. 

### Input
- directory with experimental MS files with mass and intensity columns in .txt format
- directory with theoretical MS files with mass and intensity columns in .txt format

### Output
- overview table showing number of peaks that between all experimental an theoretical MS files in .txt format
- experimental MS files with additional "subclass" column in .csv format

## Generate library files using MS-Digest
1. Find sequence of subclass on Uniprot or IMGT
2. Use [MS-Digest](https://prospector.ucsf.edu/prospector/cgi-bin/msform.cgi?form=msdigest) to create theoretical spectra
3. Collect data