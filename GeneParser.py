import os

# Create directory for processed sequences
os.makedirs('ProcessedAASequences', exist_ok=True)
foldername = 'YeastGenes'
RGCAULevels = open("RGCAULevels.txt", "w")
headers = ["Gene", "GC", "AU", "Length", "Polar%", "Non-Polar%"]
RGCAULevels.write('{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\n'.format(*headers))

# Define the codon table
codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Define amino acids for polar and non-polar
polar_amino_acids = set('DEHKNQRSTY')
non_polar_amino_acids = set('ACFGILMPVW')

def main():
    for filename in os.listdir(foldername):
        # Initialize values
        gene_name = os.path.splitext(filename)[0]
        gc_content = au_content = polar_percent = non_polar_percent = 0
        sequence = ""

        # Read sequence from file
        with open(os.path.join(foldername, filename)) as infile:
            for line in infile:
                if not line.startswith('>'):  # Skip header lines
                    sequence += line.strip()

        # Calculate GC Content
        gc_content = round((sequence.count('G') + sequence.count('C')) / len(sequence) * 100, 1)

        # Calculate AU Content
        au_content = round((sequence.count('A') + sequence.count('U')) / len(sequence) * 100, 1)

        # Calculate polar and non-polar percentages
        polar_count = sum(sequence.count(aa) for aa in polar_amino_acids)
        non_polar_count = sum(sequence.count(aa) for aa in non_polar_amino_acids)
        total_amino_acids = len(sequence)
        polar_percent = round((polar_count / total_amino_acids) * 100, 1)
        non_polar_percent = round((non_polar_count / total_amino_acids) * 100, 1)

        # Write results to file
        RGCAULevels.write('{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\t{:<10}\n'.format(gene_name, gc_content, au_content, len(sequence), polar_percent, non_polar_percent))

    RGCAULevels.close()
    print("end myGeneParser.py")

if __name__ == "__main__":
    main()