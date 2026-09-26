import sys
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

        # Keep the sequence in its original case rather than calling
        # .upper() on the whole thing. Biopython's complement() is
        # case-aware (a<->t, c<->g, case preserved), so complementing
        # doesn't require uppercasing first -- this drops one full-sequence
        # traversal that wasn't actually needed.
        #
        # NOTE: since self.seq keeps its original case, any code comparing
        # individual bases (e.g. Ursina's renderer) needs to handle lowercase
        # too, or uppercase just the small slice it's about to render.
        self.seq = record.seq
        self.complement = self.seq.complement()
        self.length = len(self.seq)

        # Count once on the original (mixed-case) sequence, then fold
        # lowercase keys into their uppercase equivalents. The fold only
        # ever touches a handful of possible characters (a/c/g/t/n and rare
        # IUPAC codes), so it's effectively free -- versus calling .upper()
        # on the whole sequence just to make counting case-insensitive,
        # which would cost a second full traversal for the same result.
        raw_counts = Counter(str(self.seq))
        self.base_counts = Counter()
        for base, count in raw_counts.items():
            self.base_counts[base.upper()] += count

        # known_bases = only the four called nucleotides, excluding N and
        # any rarer IUPAC ambiguity codes. GC% is computed over this instead
        # of `length`, so gap/unresolved regions don't quietly deflate the
        # percentage the way dividing by the full length does.
        self.known_bases = (
            self.base_counts["A"]
            + self.base_counts["C"]
            + self.base_counts["G"]
            + self.base_counts["T"]
        )
        gc = self.base_counts["G"] + self.base_counts["C"]
        self.gc_content = (gc / self.known_bases * 100) if self.known_bases else 0.0

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
        self.total_known_bases = 0
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
            self.total_known_bases += chrom.known_bases

        gc = self.total_base_counts["G"] + self.total_base_counts["C"]
        self.total_gc_content = (
            (gc / self.total_known_bases * 100) if self.total_known_bases else 0.0
        )

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

    def getTotalKnownBases(self):
        return self.total_known_bases

    def getTotalGC_Content(self):
        return self.total_gc_content

    # ---- misc ----

    def printSequenceSizesBytes(self):
        """
        Print the in-memory size, in bytes, of each chromosome's stored
        sequence string, plus a running total.

        Note: self.seq is stored as a Biopython Seq object, not a plain str,
        and sys.getsizeof() on a Seq doesn't reliably reflect the size of
        its underlying data. To measure honestly, this briefly converts each
        Seq to str with str(chrom.seq) -- a temporary, one-off copy made
        just for this measurement, not something kept around afterward. For
        very large chromosomes that temporary copy is itself a non-trivial
        chunk of memory for the moment it exists, so treat this as a
        diagnostic to run occasionally, not something to call in a loop.
        """
        total_bytes = 0
        for chrom_id in self.getChromosomeIds():
            chrom = self.getChromosome(chrom_id)
            size_bytes = sys.getsizeof(str(chrom.seq))
            total_bytes += size_bytes
            print(f"{chrom_id}: {size_bytes:,} bytes")

        print(f"Total (sum of chromosome seq strings): {total_bytes:,} bytes")

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
    print(f"Genome-wide known bases: {p.getTotalKnownBases()}")
    print(f"Genome-wide GC%: {p.getTotalGC_Content():.2f}")

    first_id = p.getChromosomeIds()[0]
    print(f"First record: {first_id}, GC% = {p.getGC_Content(first_id):.2f}")

    p.printSequenceSizesBytes()
    # p.printSequenceAndComplement("NC_000011.10")


if __name__ == "__main__":
    main()