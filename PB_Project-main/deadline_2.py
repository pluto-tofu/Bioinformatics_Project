from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW

RESOURCE_DIR = "resources"

def save_to_file(data, file_path, mode='a'):
    """Utility function to write data to a file."""
    with open(file_path, mode) as file:
        file.write(data)

def perform_blast(accession, filename, email, rice):
    """Performs a BLAST search for a given accession number."""
    try:
        Entrez.email = email
        result_handle = NCBIWWW.qblast("blastn", "nt", accession, entrez_query=f"{rice}[AC]")

        output_file = f"{RESOURCE_DIR}/{filename}.xml"
        save_to_file(result_handle.read(), output_file)
        result_handle.close()

        print(f"BLAST result saved to {output_file}")
    except Exception as e:
        print(f"An error occurred while performing BLAST: {e}")

email="rishit22405@iiitd.ac.in"
accessions = [
    ("AB021878.1", "salinity_tolerance_rice"),
    ("MH881012.1", "starch_quality_rice"),
    ("AB028184.1", "drought_resistance_rice"),
]

rice_genes={
    "AB021878.1":["SRX8723763","DRX161746"],
    "MH881012.1":["SRX20596260","SRX24433606"],
    "AB028184.1":["SRX13519067","SRX510041","SRX19861552"]
}

for accession, name in accessions:
    open(f"{RESOURCE_DIR}/{name}.xml",'w').close()
    if accession in rice_genes.keys():
        for rice in rice_genes[accession]:
            perform_blast(accession, name, email, rice)