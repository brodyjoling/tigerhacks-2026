from Bio import SeqIO
from collections import Counter
from time import perf_counter

class ProtoDNASequencer:

    def __init__(self, fasta_file):
        self.fasta_file = fasta_file
        self.record = None
        self.seq = None
        self.complement = None
        self.base_counts = None
        self.gc_content = None

    def load_sequence(self):
        start_time = perf_counter()
        self.record = next(SeqIO.parse(self.fasta_file, "fasta"))
        elapsed_time = perf_counter() - start_time
        print(f"Time to sequence FASTA file: {elapsed_time:.6f} seconds")
        self.seq = self.record.seq
        self.complement = self.seq.complement()
        self.base_counts = Counter(str(self.seq))
        self.gc_content = (self.base_counts["G"] + self.base_counts["C"]) / len(self.seq) * 100

    def getHelix(self):
        pass

    def getFile(self):
        return self.fasta_file

    def getRecord(self):
        return self.record

    def getSeq(self):
        return self.seq

    def getComplement(self):
        return self.complement

    def getBaseCounts(self):
        return self.base_counts

    def getGC_Content(self):
        return self.gc_content


def main():
    p = ProtoDNASequencer("GCF_000001405.40_GRCh38.p14_genomic.fna")
    p.load_sequence()
    print(p.getSeq())


if __name__ == "__main__":
    main()


# base_pairs = []
# for i, (base, comp_base) in enumerate(zip(str(seq), str(complement))):
#     base_pairs.append({
#         "index": i,
#         "base": base,
#         "complement": comp_base,
#     })

# def printPairsAndInfo(pairs):
#     for x in pairs:
#         print(str(x.get("index")) + ": " + x.get("base") + x.get("complement"))
#         # sleep(0.001)

# printPairsAndInfo(base_pairs)