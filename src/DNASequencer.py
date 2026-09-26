from pyfaidx import Fasta
from collections import Counter
import numpy as np
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
fasta_path = os.path.join(script_dir, "GCF_000001405.40_GRCh38.p14_genomic.fna")

genome = Fasta(fasta_path)

# See what's in the file without loading any sequence data
print("Chromosomes/records found:")
for name in genome.keys():
    print(f"  {name}: {len(genome[name])} bases")

def get_window(chrom_name, start, length=200):
    """Fetch a small slice of sequence without loading the whole chromosome."""
    region = genome[chrom_name][start:start+length]
    return str(region)

def complement_str(seq_str):
    pairing = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(pairing.get(b, 'N') for b in seq_str)

def get_base_pairs(chrom_name, start, length=200):
    """Build the per-base dict structure, but only for a small window."""
    window = get_window(chrom_name, start, length)
    comp = complement_str(window)
    return [
        {"index": start + i, "base": b, "complement": c}
        for i, (b, c) in enumerate(zip(window, comp))
    ]

# Example: grab a 200-base window starting at position 5,225,464 on chromosome 11
pairs = get_base_pairs("NC_000011.10", 5225464, 200)
for p in pairs[:10]:
    print(f"{p['index']}: {p['base']}{p['complement']}")

def gc_content(seq_str):
    counts = Counter(seq_str)
    total = len(seq_str)
    return (counts.get("G", 0) + counts.get("C", 0)) / total * 100 if total else 0

window = get_window("NC_000011.10", 5225464, 200)
print(f"GC content of this window: {gc_content(window):.1f}%")