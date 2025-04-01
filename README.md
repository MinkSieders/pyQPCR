# pyQPCR

## Description
A script for analyzing qPCR results, processing raw spectra, dissociation curves, and performing qPCR plate design. It can analyze Ct values, spectra data, and dissociation curves for each well of a qPCR plate.

## Main Commands
**Analyzer** (_analyzer_): Analyzes qPCR results from input files containing Ct, spectra, and dissociation curve data.

**Plater** (_plater_): Designs a qPCR plate based on the sample list or layout file.

## Input
Each command is accessed through a command-line argument specifying the mode, followed by additional flags for customization. Most commands are fitted with default values; make sure to change according to your needs.

For in-depth explanation on all the flags within each command, use --help or -h.

## Usage | Plater

`python pyQPCR.py plater`

The plater command allows users to design a qPCR plate layout based on a sample list or prepared layout file.

### Arguments:

`--input` (str, required): Path to the sample list (.txt) or prepared layout file (.tsv).

`--design` (str, optional, default: sequential_1): Plate design for automatic well assignment. Currently supports sequential_1.

`--target_genes` (list of str, optional): List containing 3-letter primer gene target IDs.

`--number_technical_replicates` (int, optional, default: 2): Number of technical replicates for each sample.

### Example Usage:

To design a qPCR plate with a sequential layout of 96 wells and a custom list of target genes:

`python pyQPCR.py plater --input samples.txt --design sequential_1 --target_genes GeneX GeneY GeneZ --number_technical_replicates 3`

## Output Handling

For each command different outputs are generated based on each different analysis. Output file paths are determined based in the input file path name with extensions based on the type of analysis / processing performed. 

## Author
Mink Sieders

## License
MIT License. Copyright (c) 2024 Mink Sieders

