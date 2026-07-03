#!/bin/bash
echo "Starting DNA Sequence Analysis..."
echo "Input file: data/sequences.fasta"
python3 analyze.py --input data/sequences.fasta
echo "Analysis complete! Check outputs/ for results."
