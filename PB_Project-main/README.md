# PB_Project

## Description

This Python script automates the downloading of nucleotide sequences from the NCBI database and performs BLAST searches against the "nt" database, specifically filtering results by the organism "Oryza sativa". The script is designed for researchers needing automated, batch downloads of sequences and their subsequent analysis.

We performed test on the rice genome and on different strains of the characteristics that different varieties of rice show
out prime focus is on how to pics the perfect rice which can give us better starch quality and shows salt tolerance or drought tolerance as needed for a particular area


## Features

- **Download Sequences:** Fetches FASTA formatted sequences from NCBI.
- **BLAST Search:** Performs a BLAST search on downloaded sequences.
- **Save Results:** Saves both the FASTA files and BLAST results locally for further analysis.

## Usage


### Running the Script

To run the script, use the following command from the terminal:

```bash
python script_name.py
```
### Output

The script will download FASTA files and BLAST XML results to the `resources` directory. Each file is named according to the sequence's accession number and the type of analysis, e.g., `salinity_tolerance.fasta` or `salinity_tolerance.xml`.

## Customization

You can modify the script to add more sequences or change the target organism by editing the `sequences` tuple in the `main` function:

```python
sequences = [
    ("New_Accession_Number", "nuccore", "Descriptive_Name"),
    ...
]
```
## Troubleshooting

Ensure that your internet connection is stable during the operation as the script requires access to NCBI's servers. If you encounter any errors related to NCBI data fetching, verify that your Entrez email is correctly set and that the accession numbers are valid.


## References

### Starch Quality:
Modification of Waxy Gene improves Starch Quality.

- **Research Paper:** https://pubmed.ncbi.nlm.nih.gov/34063649/

- **mRNA:** https://www.ncbi.nlm.nih.gov/nuccore/MH881012.1?report=fasta

### Salt Tolerance:
OsNHX1 gene is responsible for Salt Tolerance.

- **Research Paper:** https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4724728/OsNHX1

- **mRNA:** https://www.ncbi.nlm.nih.gov/nuccore/AB021878.1?report=fasta

### Drought Resistance:
Absence of synaptotagmin-5 (OsSYT-5) highly increase the Drought Resistance.

- **Research Paper:** https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0258171

- **mRNA:** https://www.ncbi.nlm.nih.gov/nuccore/XM_015795353.2?report=fasta

- **Protein:** https://www.ncbi.nlm.nih.gov/protein/XP_025878309.1?report=fasta

## We have performed BLAST on these FASTA files to determine the DNA regions where these genes reside.
