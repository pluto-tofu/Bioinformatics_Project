from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW

# Constants
RESOURCE_DIR = "resources"

def save_to_file(data, file_path, mode='w'):
    """Utility function to write data to a file."""
    with open(file_path, mode) as file:
        file.write(data)

def download_fasta(accession, db, filename, email):
    """Downloads a FASTA file from NCBI."""
    try:
        Entrez.email = email
        handle = Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()

        output_file = f"{RESOURCE_DIR}/{filename}.fasta"
        with open(output_file, "w") as file:
            SeqIO.write(record, file, "fasta")

        print(f"FASTA sequence saved to {output_file}")
    except Exception as e:
        print(f"An error occurred while downloading FASTA: {e}")

def perform_blast(accession, filename, email):
    """Performs a BLAST search for a given accession number."""
    try:
        Entrez.email = email
        result_handle = NCBIWWW.qblast("blastn", "nt", accession, entrez_query="Oryza sativa[organism]")

        output_file = f"{RESOURCE_DIR}/{filename}.xml"
        save_to_file(result_handle.read(), output_file)
        result_handle.close()

        print(f"BLAST result saved to {output_file}")
    except Exception as e:
        print(f"An error occurred while performing BLAST: {e}")

# Setup parameters
email = "satwik22461@iiitd.ac.in"
accessions = [
    ("AB021878.1", "salinity_tolerance"),
    ("MH881012.1", "starch_quality"),
    ("AB028184.1", "drought_resistance"),
]

# Download FASTA sequences
for accession, name in accessions:
    download_fasta(accession, "nuccore", name, email)

# Perform BLAST analyses
for accession, name in accessions:
    perform_blast(accession, name, email)

