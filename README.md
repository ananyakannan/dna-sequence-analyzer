# DNA Sequence Analyzer
A python bioinformatics tool that analyzes DNA sequences from FASTA files. Built using Biopython, Pandas and Matplotlib.

## What it does
- Parses FASTA file using Biopython
- Calcutes GC content, AT content, and GC/AT ratio
- Counts nucleotides composition (A, T, G, C)
- Detects Open Reading Frames (ORFs)
- Exports resutls to CSV
- Generates 4 publication-style plots
- Flags sequences with no ORFs as a warning

# How to run
```bash
bash run_analysis.sh
```

Or directly with Python:

```bash
python3 analyze.py --input data/sequences.fasta
```
## Output
- `outputs/results.csv` — full analysis table
- `outputs/plots/gc_content_bar.png` — GC content bar chart
- `outputs/plots/length_histogram.png` — sequence length distribution
- `outputs/plots/nucleotide_composition.png` — nucleotide stacked bar chart
- `outputs/plots/gc_vs_length_scatter.png` — GC% vs length scatter plot
- `outputs/plots/combined_plots.png` — summary panel with all 4 plots

## Tools used

- Python 3
- Biopython
- Pandas
- Matplotlib

## Author

Ananya Kannan - MS Biotechnology, Northeastern University
