#analyze .py
#DNA sequence analyzer- Day 2
from Bio import SeqIO
fasta_file = "data/sequences.fasta"
def calculate_gc_content(sequence):
    sequence = str(sequence)
    gc = (sequence.count("G") + sequence.count("C")) / len(sequence) * 100
    return round(gc, 2)
def count_nucleotides(sequence):
    sequence = str(sequence)
    counts = {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C")
    }
    return counts
def get_sequence_length(sequence):
    return len(sequence)
for record in SeqIO.parse(fasta_file, "fasta"):
    gc = calculate_gc_content(record.seq)
    counts = count_nucleotides(record.seq)
    length = get_sequence_length(record.seq)
    print(f" ID: {record.id}")
    print(f" Length: {length}")
    print(f" GC Content: {gc}%")
    print(f" Nucleotide Counts: {counts}")
    print("---")



    

