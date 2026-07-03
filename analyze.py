#analyze.py
#DNA sequence Analyzer

from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="DNA Sequence Analyzer")
parser.add_argument("--input", required=True, help="Path to input FASTA file")
args = parser.parse_args()
fasta_file = args.input

def calculate_gc_content(sequence):
    sequence = str(sequence)
    gc = (sequence.count("G") + sequence.count("C")) / len(sequence) * 100
    return round(gc,2)

def count_nucleotides(sequence):
    sequence = str(sequence)
    counts = {
        "A" : sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C")
    }
    return counts

def get_sequence_lengths(sequence):
    return len(sequence)

def find_orfs(sequence):
    sequence = str(sequence)
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []
    for i in range(len(sequence)):
        if sequence[i:i+3] == start_codon:
            for j in range(i, len(sequence), 3):
                codon = sequence[j:j+3]
                if codon in stop_codons:
                    orfs.append(sequence[i:j+3])
                    break
    return orfs

results = []
for record in SeqIO.parse(fasta_file, "fasta"):
    gc = calculate_gc_content(record.seq)
    counts = count_nucleotides(record.seq)
    length = get_sequence_lengths(record.seq)
    orfs = find_orfs(record.seq)
    row = {
        "id": record.id,
        "length": length,
        "gc_content": gc,
        "A": counts["A"],
        "T": counts["T"],
        "G": counts["G"],
        "C": counts["C"], 
        "orf_counts": len(orfs)
    }
    results.append(row)

df = pd.DataFrame(results)
df["at_content"] = round(100-df["gc_content"], 2)
df["gc_at_ratio"] = round(df["gc_content"]/df["at_content"], 2)
no_orf_sequences = df[df["orf_counts"] == 0]
if len(no_orf_sequences) > 0:
    print("WARNING: These sequences have no ORFs:")
    print(no_orf_sequences["id"].values)
df.to_csv("outputs/results.csv", index=False)
print("Results saved to outputs/results.csv")

plt.figure(figsize=(10, 6))
plt.bar(df["id"], df["gc_content"], color="steelblue")
plt.title("GC Content per Sequence")
plt.xlabel("Sequence ID")
plt.ylabel("GC Content (%)")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("outputs/plots/gc_content_bar.png")
plt.close()
print("GC content bar chart saved!")

plt.figure(figsize=(10, 6))
plt.hist(df["length"], bins=5, color="green", edgecolor="black")
plt.title("Sequence Length Distribution")
plt.xlabel("Sequence Length (bp)")
plt.ylabel("Number of Sequences")
plt.tight_layout()
plt.savefig("outputs/plots/length_histogram.png")
plt.close()
print("Length histogram saved!")

plt.figure(figsize=(10, 6))
plt.bar(df["id"], df["A"], label="A", color="steelblue")
plt.bar(df["id"], df["T"], bottom=df["A"], label="T", color="coral")
plt.bar(df["id"], df["G"], bottom=df["A"]+df["T"], label="G", color="mediumseagreen")
plt.bar(df["id"], df["C"], bottom=df["A"]+df["T"]+df["G"], label="C", color="gold")
plt.title("Nucleotide Composition per Sequence")
plt.xlabel("Sequence ID")
plt.ylabel("Nucleotide Count")
plt.legend()
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("outputs/plots/nucleotide_composition.png")
plt.close()
print("Nucleotide composition chart saved!")

plt.figure(figsize=(10, 6))
plt.scatter(df["length"], df["gc_content"], color="purple", s=100, zorder=5)
plt.axhline(y=50, color="red", linestyle="--", label="50% GC reference")
plt.title("GC Content vs Sequence Length")
plt.xlabel("Sequence Length (bp)")
plt.ylabel("GC Content (%)")
plt.legend()
plt.tight_layout()
plt.savefig("outputs/plots/gc_vs_length_scatter.png")
plt.close()
print("Scatter plot saved!")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes[0][0].bar(df["id"], df["gc_content"], color="steelblue")
axes[0][0].set_title("GC Content per Sequence")
axes[0][0].tick_params(axis='x', rotation=45)
axes[0][1].hist(df["length"], bins=5, color="green", edgecolor="black")
axes[0][1].set_title("Sequence Length Distribution")
axes[1][0].bar(df["id"], df["A"], label="A", color="steelblue")
axes[1][0].bar(df["id"], df["T"], bottom=df["A"], label="T", color="coral")
axes[1][0].bar(df["id"], df["G"], bottom=df["A"]+df["T"], label="G", color="mediumseagreen")
axes[1][0].bar(df["id"], df["C"], bottom=df["A"]+df["T"]+df["G"], label="C", color="gold")
axes[1][0].set_title("Nucleotide Composition")
axes[1][0].tick_params(axis='x', rotation=45)
axes[1][0].legend()
axes[1][1].scatter(df["length"], df["gc_content"], color="purple", s=100)
axes[1][1].axhline(y=50, color="red", linestyle="--")
axes[1][1].set_title("GC Content vs Sequence Length")
plt.tight_layout()
plt.savefig("outputs/plots/combined_plots.png")
plt.close()
print("Summary panel saved!")

print("\n=== ANALYSIS SUMMARY ===")
print(f"Total sequences analyzed: {len(df)}")
print(f"Average GC content: {round(df['gc_content'].mean(), 2)}%")
print(f"Average sequence length: {round(df['length'].mean(), 2)}bp")
print(f"Total ORFs found: {df['orf_counts'].sum()}")
print(f"Sequences with no ORFs: {len(df[df['orf_counts'] == 0])}")
print("=============================")







    

