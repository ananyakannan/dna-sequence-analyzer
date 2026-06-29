#analyze .py
#DNA sequence analyzer- Day 1
from Bio import SeqIO
#This is the path to our FASTA file
fasta_file = "data/sequences.fasta"
#SeqIO.parse() reads the file and gives us one Seqeunce Record at a time
#"fasta" tells in what format to expect
for record in SeqIO.parse(fasta_file, "fasta"):
    length = len(record.seq)
    print(f"ID: {record.id} | Length: {length} bp")

