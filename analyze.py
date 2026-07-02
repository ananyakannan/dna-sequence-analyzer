#analyze .py
#DNA sequence analyzer- Day 2
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
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
results = []
for record in SeqIO.parse(fasta_file, "fasta"):
    gc = calculate_gc_content(record.seq)
    counts = count_nucleotides(record.seq)
    length = get_sequence_length(record.seq)
    row = {
    "id": record.id,
    "length": length,
    "gc_content": gc,
    "A": counts["A"],
    "T": counts["T"],
    "C": counts["C"],
    "G": counts["G"]
}
    results.append(row)
df = pd.DataFrame(results)
df["at_content"] = round(100 - df["gc_content"], 2)
df["gc_at_ratio"] = round(df["gc_content"] / df["at_content"], 2)
df.to_csv("outputs/results.csv", index=False)
print("Results saved to outputs/results.csv")
plt.figure(figsize=(10, 6))
plt.bar(df["id"], df["gc_content"], color="steelblue")
plt.title("GC Content per Sequence")
plt.xlabel("Sequence ID")
plt.ylabel("GC Content (%)")
plt.xticks(rotation=45)
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
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes[0][0].bar(df["id"], df["gc_content"], color="steelblue")
axes[0][0].set_title("GC content per Sequence")
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
print("GC content vs sequence length scatter plot saved!")










    

