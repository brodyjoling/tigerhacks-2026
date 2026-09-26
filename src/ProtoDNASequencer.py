from Bio import SeqIO
from collections import Counter
from time import perf_counter


class ChromosomeData:
    """
    Holds everything we know about ONE chromosome/contig record from the
    FASTA file (its own seq, complement, base_counts, gc_content, etc.).

    Rationale: rather than keeping four parallel lists on the main sequencer
    (all_seqs, all_complements, all_base_counts, all_gc_contents) that have
    to stay in sync by index, each chromosome's data lives together in one
    small object. "Give me chromosome 7's info" becomes a single lookup
    instead of hunting through several separate collections for the same
    index.
    """

    def __init__(self, record):
        self.record = record
        self.id = record.id                # e.g. "NC_000001.11"
        self.description = record.description

        # Normalize case once, up front. Some assemblies soft-mask repeat
        # regions in lowercase (a/c/g/t); without this, a case-sensitive
        # Counter undercounts G/C in those regions while len() still counts
        # them, silently deflating gc_content. Uppercasing here means every
        # attribute below (complement, base_counts, gc_content) is
        # automatically consistent -- no case-insensitive counting logic
        # needed anywhere else.
        self.seq = record.seq.upper()
        self.complement = self.seq.complement()
        self.base_counts = Counter(str(self.seq))
        self.length = len(self.seq)
        gc = self.base_counts["G"] + self.base_counts["C"]
        self.gc_content = (gc / self.length * 100) if self.length else 0.0

    def __repr__(self):
        return f"<ChromosomeData {self.id} len={self.length} GC={self.gc_content:.2f}%>"


class ProtoDNASequencer:
    """
    Parses every record (chromosome/contig) in a FASTA file.

    Storage strategy
    -----------------
    self.chromosomes : dict[str, ChromosomeData]
        Keyed by record.id. Lets you pull a single chromosome's
        seq/complement/base_counts/gc_content directly, e.g.
        sequencer.chromosomes["NC_000001.11"].gc_content

    self.total_base_counts / self.total_length / self.total_gc_content
        Genome-wide numbers, built by ACCUMULATING each chromosome's counts
        as we go, rather than concatenating every chromosome's sequence and
        complement into one giant string and recomputing from that.

    Why accumulate instead of concatenate
    --------------------------------------
    GRCh38 is ~3.1 billion bases. A literal genome-wide `seq` string plus
    its `complement` string, on top of the per-chromosome copies already
    held in `self.chromosomes`, means multiple full copies of the genome
    in memory at once. Base composition and GC% are additive
    (total_G = sum of each chromosome's G count, etc.), so a running
    Counter gives the exact same genome-wide statistics for a small,
    constant amount of extra memory -- no second copy of the genome
    required.

    If a literal, continuous, concatenated genome sequence turns out to be
    needed later (as opposed to just genome-wide stats), that's a separate,
    heavier design question -- worth revisiting deliberately rather than
    defaulting into it here.
    """

    def __init__(self, fasta_file):
        self.fasta_file = fasta_file
        self.chromosomes = {}          # record.id -> ChromosomeData
        self.total_base_counts = Counter()
        self.total_length = 0
        self.total_gc_content = 0.0

    def load_sequence(self, chromosomes_only=True):
        """
        chromosomes_only=True (default): only keep records whose accession
        starts with "NC_" -- RefSeq's prefix for complete chromosome-level
        molecules (22 autosomes + X + Y + the mitochondrial genome, 25
        records for GRCh38.p14). Everything else in the file (~680 records,
        prefixed NT_/NW_) is unplaced/unlocalized scaffolds and patch
        sequences -- skipped here rather than stored and then ignored later.

        Pass chromosomes_only=False to keep every record, including
        scaffolds/patches, if you end up needing them.
        """
        start_time = perf_counter()

        for record in SeqIO.parse(self.fasta_file, "fasta"):
            if chromosomes_only and not record.id.startswith("NC_"):
                continue

            chrom = ChromosomeData(record)
            self.chromosomes[chrom.id] = chrom

            # Roll this chromosome's numbers into the genome-wide totals.
            self.total_base_counts.update(chrom.base_counts)
            self.total_length += chrom.length

        gc = self.total_base_counts["G"] + self.total_base_counts["C"]
        self.total_gc_content = (gc / self.total_length * 100) if self.total_length else 0.0

        elapsed_time = perf_counter() - start_time
        print(f"Time to parse {len(self.chromosomes)} sequence(s): {elapsed_time:.6f} seconds")

    def getHelix(self):
        pass

    def getFile(self):
        return self.fasta_file

    # ---- per-chromosome accessors ----

    def getChromosomeIds(self):
        """List every chromosome/contig id that was parsed."""
        return list(self.chromosomes.keys())

    def getChromosome(self, chrom_id):
        """Return the ChromosomeData for one chromosome/contig by its record id."""
        return self.chromosomes.get(chrom_id)

    def getRecord(self, chrom_id):
        chrom = self.getChromosome(chrom_id)
        return chrom.record if chrom else None

    def getSeq(self, chrom_id):
        chrom = self.getChromosome(chrom_id)
        return chrom.seq if chrom else None

    def getComplement(self, chrom_id):
        chrom = self.getChromosome(chrom_id)
        return chrom.complement if chrom else None

    def getBaseCounts(self, chrom_id):
        chrom = self.getChromosome(chrom_id)
        return chrom.base_counts if chrom else None

    def getGC_Content(self, chrom_id):
        chrom = self.getChromosome(chrom_id)
        return chrom.gc_content if chrom else None

    # ---- whole-genome accessors ----

    def getTotalBaseCounts(self):
        return self.total_base_counts

    def getTotalLength(self):
        return self.total_length

    def getTotalGC_Content(self):
        return self.total_gc_content

    # ---- misc ----

    def printSequenceAndComplement(self, chrom_id):
        chrom = self.getChromosome(chrom_id)
        if chrom is None:
            print(f"No chromosome loaded with id '{chrom_id}'.")
            return

        print(f"Index\tBase\tComplement  ({chrom_id})")
        for i, (base, comp_base) in enumerate(zip(str(chrom.seq), str(chrom.complement))):
            print(f"{i}\t{base}\t{comp_base}")


def main():
    p = ProtoDNASequencer("data/GCF_000001405.40_GRCh38.p14_genomic.fna")
    p.load_sequence()

    print(f"Parsed {len(p.getChromosomeIds())} sequence(s)")
    print(f"Genome-wide length: {p.getTotalLength()}")
    print(f"Genome-wide GC%: {p.getTotalGC_Content():.2f}")

    first_id = p.getChromosomeIds()[0]
    print(f"First record: {first_id}, GC% = {p.getGC_Content(first_id):.2f}")


if __name__ == "__main__":
    main()