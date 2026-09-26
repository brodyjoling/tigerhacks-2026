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

        # --- GC / base-count logic disabled for now ---
        # Not needed yet, and skipping it avoids a full traversal of the
        # (potentially huge) sequence every time a chromosome is sequenced.
        # Re-enable if/when GC% or base composition is actually needed.
        #
        # # Count once on the original (mixed-case) sequence, then fold
        # # lowercase keys into their uppercase equivalents. The fold only
        # # ever touches a handful of possible characters (a/c/g/t/n and rare
        # # IUPAC codes), so it's effectively free -- versus calling .upper()
        # # on the whole sequence just to make counting case-insensitive,
        # # which would cost a second full traversal for the same result.
        # raw_counts = Counter(str(self.seq))
        # self.base_counts = Counter()
        # for base, count in raw_counts.items():
        #     self.base_counts[base.upper()] += count
        #
        # # known_bases = only the four called nucleotides, excluding N and
        # # any rarer IUPAC ambiguity codes. GC% is computed over this instead
        # # of `length`, so gap/unresolved regions don't quietly deflate the
        # # percentage the way dividing by the full length does.
        # self.known_bases = (
        #     self.base_counts["A"]
        #     + self.base_counts["C"]
        #     + self.base_counts["G"]
        #     + self.base_counts["T"]
        # )
        # gc = self.base_counts["G"] + self.base_counts["C"]
        # self.gc_content = (gc / self.known_bases * 100) if self.known_bases else 0.0

    def __repr__(self):
        return f"<ChromosomeData {self.id} len={self.length}>"


class ProtoDNASequencer:
    """
    Indexes every record (chromosome/contig) in a FASTA file, but does NOT
    read any sequence data into memory until sequenceChromosome() is told
    which chromosome to actually sequence.

    Storage strategy
    -----------------
    self._index : Bio.File._IndexedSeqFileDict (from SeqIO.index)
        A lazy, dict-like index of record id -> file offset. Building this
        only reads headers/offsets, not sequence data, so it stays cheap
        even for a ~3.1 billion base genome.

    self.chromosome_ids : list[str]
        Every id available to sequence (after chromosomes_only filtering),
        whether or not it has actually been sequenced yet.

    self.chromosomes : dict[str, ChromosomeData]
        Keyed by record.id, but only populated for chromosomes that have
        actually been passed to sequenceChromosome(). Lets you pull a single
        sequenced chromosome's seq/complement directly, e.g.
        sequencer.chromosomes["NC_000001.11"].seq

    Why lazy instead of parsing everything up front
    --------------------------------------------------
    GRCh38 is ~3.1 billion bases. Eagerly parsing every chromosome means
    paying that memory/CPU cost even if the caller only ever wants one or
    two chromosomes. Indexing first and sequencing on demand means the cost
    of a chromosome is only paid when something actually asks for it.

    NOTE: total_base_counts / total_gc_content logic has been disabled for
    now (see comments below) -- it doesn't make much sense to track
    genome-wide totals when chromosomes are sequenced one at a time anyway.
    """

    def __init__(self, fasta_file):
        self.fasta_file = fasta_file
        self._index = None              # SeqIO.index(...) -- built by load_sequence()
        self.chromosome_ids = []        # ids available to sequence (post-filtering)
        self.chromosomes = {}           # record.id -> ChromosomeData, populated lazily

        self.total_length = 0           # running total over chromosomes actually sequenced so far

        # --- genome-wide base-count / GC logic disabled for now ---
        # self.total_base_counts = Counter()
        # self.total_known_bases = 0
        # self.total_gc_content = 0.0

    def load_sequence(self, chromosomes_only=True):
        """
        Index every record in the FASTA file -- id and file offset only, no
        sequence data read yet -- then narrow self.chromosome_ids down to
        the ones this sequencer will allow sequencing.

        chromosomes_only=True (default): only keep records whose accession
        starts with "NC_" -- RefSeq's prefix for complete chromosome-level
        molecules (22 autosomes + X + Y + the mitochondrial genome, 25
        records for GRCh38.p14). Everything else in the file (~680 records,
        prefixed NT_/NW_) is unplaced/unlocalized scaffolds and patch
        sequences -- excluded here rather than made sequenceable and then
        ignored later.

        Pass chromosomes_only=False to allow sequencing any record,
        including scaffolds/patches, if you end up needing them.

        Call sequenceChromosome(chrom_id) afterward to actually parse and
        store one chromosome at a time.
        """
        start_time = perf_counter()

        self._index = SeqIO.index(self.fasta_file, "fasta")

        if chromosomes_only:
            self.chromosome_ids = [rid for rid in self._index.keys() if rid.startswith("NC_")]
        else:
            self.chromosome_ids = list(self._index.keys())

        elapsed_time = perf_counter() - start_time
        print(f"Time to index {len(self.chromosome_ids)} sequence(s): {elapsed_time:.6f} seconds")

    def sequenceChromosome(self, chrom_id):
        """
        Actually parse ONE chromosome and store its sequence + complement.

        This is the only method that reads sequence data into memory. Every
        other accessor just reads from whatever has already been sequenced.
        If chrom_id was already sequenced, returns the cached ChromosomeData
        instead of re-parsing it.
        """
        if chrom_id in self.chromosomes:
            return self.chromosomes[chrom_id]

        if self._index is None:
            raise RuntimeError("Call load_sequence() before sequencing a chromosome.")

        if chrom_id not in self.chromosome_ids:
            if chrom_id in self._index:
                print(f"'{chrom_id}' exists in the file but was excluded by "
                      f"chromosomes_only filtering in load_sequence().")
            else:
                print(f"No record with id '{chrom_id}' found in {self.fasta_file}.")
            return None

        start_time = perf_counter()

        record = self._index[chrom_id]       # only this record's data is read here
        chrom = ChromosomeData(record)
        self.chromosomes[chrom_id] = chrom
        self.total_length += chrom.length

        # --- genome-wide base-count / GC accumulation disabled for now ---
        # self.total_base_counts.update(chrom.base_counts)
        # self.total_known_bases += chrom.known_bases
        # gc = self.total_base_counts["G"] + self.total_base_counts["C"]
        # self.total_gc_content = (
        #     (gc / self.total_known_bases * 100) if self.total_known_bases else 0.0
        # )

        elapsed_time = perf_counter() - start_time
        print(f"Time to sequence '{chrom_id}': {elapsed_time:.6f} seconds")

        return chrom

    def sequenceALL(self, chromosomes):
        for x in chromosomes:
            self.sequenceChromosome(x)

    def getHelix(self):
        pass

    def getFile(self):
        return self.fasta_file

    # ---- per-chromosome accessors (read from cache only -- do not sequence) ----

    def getChromosomeIds(self):
        """Every chromosome id available to sequence, per load_sequence()'s filtering."""
        return list(self.chromosome_ids)

    def getSequencedChromosomeIds(self):
        """Chromosome ids that have actually been sequenced (parsed) so far."""
        return list(self.chromosomes.keys())

    def getChromosome(self, chrom_id):
        """Return the ChromosomeData for a chromosome that's already been sequenced, or None."""
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

    # ---- base-count / GC accessors disabled for now ----

    # def getBaseCounts(self, chrom_id):
    #     chrom = self.getChromosome(chrom_id)
    #     return chrom.base_counts if chrom else None

    # def getGC_Content(self, chrom_id):
    #     chrom = self.getChromosome(chrom_id)
    #     return chrom.gc_content if chrom else None

    # ---- whole-genome accessors ----

    def getTotalLength(self):
        """Total length over chromosomes sequenced so far (not the whole genome unless all have been sequenced)."""
        return self.total_length

    # def getTotalBaseCounts(self):
    #     return self.total_base_counts

    # def getTotalKnownBases(self):
    #     return self.total_known_bases

    # def getTotalGC_Content(self):
    #     return self.total_gc_content

    # ---- misc ----

    def printSequenceSizesBytes(self):
        """
        Print the in-memory size, in bytes, of each already-sequenced
        chromosome's stored sequence string, plus a running total.

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
        for chrom_id in self.getSequencedChromosomeIds():
            chrom = self.getChromosome(chrom_id)
            size_bytes = sys.getsizeof(str(chrom.seq))
            total_bytes += size_bytes
            print(f"{chrom_id}: {size_bytes:,} bytes")

        print(f"Total (sum of sequenced chromosome seq strings): {total_bytes:,} bytes")

    def printSequenceAndComplement(self, chrom_id):
        # Goes through sequenceChromosome() rather than getChromosome() so
        # this works even if chrom_id hasn't been sequenced yet -- it's the
        # one print/debug helper that actually needs the data to exist.
        chrom = self.sequenceChromosome(chrom_id)
        if chrom is None:
            return

        print(f"Index\tBase\tComplement  ({chrom_id})")
        for i, (base, comp_base) in enumerate(zip(str(chrom.seq), str(chrom.complement))):
            print(f"{i}\t{base}\t{comp_base}")


def main():
    p = ProtoDNASequencer("data/GCF_000001405.40_GRCh38.p14_genomic.fna")
    p.load_sequence()

    print(f"Indexed {len(p.getChromosomeIds())} sequence(s)")
    print(f"Sequenced so far: {len(p.getSequencedChromosomeIds())} sequence(s)")

    first_id = p.getChromosomeIds()[0]
    chrom = p.sequenceChromosome(first_id)
    print(f"Sequenced: {chrom}")

    print(f"Total length sequenced so far: {p.getTotalLength()}")

    p.printSequenceSizesBytes()
    p.printSequenceAndComplement(first_id)


if __name__ == "__main__":
    main()