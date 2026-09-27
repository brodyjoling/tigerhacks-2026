"""
Codon-level mutation comparison.

Answers: "at a given codon location (from an annotated GTF), does our
parsed reference sequence match a known ClinVar mutation, match the
healthy reference, or match neither?"

Pipeline
--------
1. ClinVar reader        -> MutationRecord (gene, transcript, genomic
                             position, ref/alt allele, protein change).
2. AnnotationLookup       -> ABSTRACT interface. Given a transcript + codon
   (ABC)                    number, returns the codon's genomic coordinates
                             and strand.
     - PostgresAnnotationLookup: the real implementation, backed by your
       teammate's Postgres 'annotations' table (find_transcripts/find_exons).
     - MockAnnotationLookup: tiny in-memory stand-in, kept around for
       offline dev/tests that shouldn't need a live DB connection.
3. CodonComparer          -> pulls the 3 bases at those coordinates out of
                             a parsed ChromosomeData (from ProtoDNASequencer),
                             applies ClinVar's ref->alt substitution, and
                             reports HEALTHY / MUTATED / UNEXPECTED.

What ClinVar actually gives you
--------------------------------
Not a literal "mutated codon" field. You get:
  - VCF (clinvar.vcf.gz): genomic position + REF/ALT allele (always
    relative to the genome's + strand, regardless of gene strand), plus
    INFO fields like CLNSIG (significance), CLNDN (condition), GENEINFO.
  - variant_summary.txt.gz: the same info in a flat TSV, plus a `Name`
    column carrying HGVS notation, e.g. "NM_000518.5(HBB):c.20A>T
    (p.Glu7Val)" -- gene, transcript, coding change, and amino acid change
    all in one string when ClinVar has classified the protein effect.
The mutant codon itself is derived, not given: reference codon (from your
parsed sequence) + ClinVar's ref->alt substitution = mutant codon.

Worked example throughout: the sickle-cell mutation.
HBB (hemoglobin beta), NM_000518.5, c.20A>T, p.Glu7Val, rs334.
HBB is on the MINUS strand of chromosome 11 (GRCh38 g.5227002T>A). Note the
flip: genomic +strand is T>A, transcript/coding strand is A>T -- exactly
the complement relationship ChromosomeData.complement already captures.

NOTE ON FIELD NAMES: the VCF INFO field names and TSV column names below
(GENEINFO, CLNSIG, CLNDN, CLNHGVS, GeneSymbol, PositionVCF, etc.) reflect
ClinVar's schema as commonly documented -- ClinVar has changed its schema
before, so cross-check against the README(_VCF).txt shipped alongside
whatever release you actually download before trusting field names blindly.

Install: pip install biopython psycopg
(pandas is no longer required -- read_clinvar_variant_summary streams the
TSV with the csv module instead; see its docstring for why.)
"""

from __future__ import annotations

import sys
from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional

from Bio.Seq import Seq

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CLINVAR_VARIANT_SUMMARY = "data/variant_summary.txt.gz"
CLINVAR_VCF = "data/clinvar.vcf.gz"
sys.path.append(str(PROJECT_ROOT))

from src.ProtoDNASequencer import ProtoDNASequencer

# Adjust this import to match wherever your teammate's DB module actually
# lives -- this assumes a module exposing find_transcripts()/find_exons()
# (the same file with find_gene_by_name/find_gene/etc.).
from parser import find_transcripts, find_exons, find_cds


# ---------------------------------------------------------------------------
# 1. ClinVar side: parsed mutation records
# ---------------------------------------------------------------------------

@dataclass
class MutationRecord:
    """One ClinVar variant, normalized to the fields this module needs."""
    gene: str
    transcript_id: str                    # e.g. "NM_000518.5"
    chrom: str                            # e.g. "11" -- normalize to match
                                           # however ProtoDNASequencer keys
                                           # its chromosomes (NC_... vs "11")
    genomic_position: int                 # 1-based, VCF convention
    ref_allele: str                       # always genome + strand
    alt_allele: str                       # always genome + strand
    protein_change: Optional[str] = None  # e.g. "p.Glu7Val", when present
    clinical_significance: Optional[str] = None
    condition: Optional[str] = None
    assembly: Optional[str] = None        # "GRCh37" or "GRCh38" (variant_summary.txt only)
    review_status: Optional[str] = None   # ClinVar's confidence tier, e.g.
                                           # "reviewed by expert panel" (variant_summary.txt only)


# GRCh38.p14 primary-assembly RefSeq accessions, keyed by the bare
# chromosome number/letter ClinVar's variant_summary.txt uses in its
# "Chromosome" column ("11", "X", "MT", etc.) -- maps to the accession
# ProtoDNASequencer actually uses as a chromosome id ("NC_000011.10").
# Verified against NCBI's GRCh38 assembly records.
GRCH38_CHROM_ACCESSIONS = {
    "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12",
    "4": "NC_000004.12", "5": "NC_000005.10", "6": "NC_000006.12",
    "7": "NC_000007.14", "8": "NC_000008.11", "9": "NC_000009.12",
    "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10",
    "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10",
    "19": "NC_000019.10", "20": "NC_000020.11", "21": "NC_000021.9",
    "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10",
    "MT": "NC_012920.1",
}


def read_clinvar_vcf(vcf_path: str):
    """
    Parse a ClinVar VCF (clinvar.vcf.gz) and yield MutationRecord objects.

    Pure-Python line parser -- no compiled dependency. cyvcf2/pysam are
    faster but both wrap htslib in C, which needs a full C build toolchain
    to install from source on Windows (no prebuilt Windows wheels exist).
    VCF is plain tab-delimited text, so a manual parser is genuinely fine
    here; it'll just be slower than cyvcf2 on the full multi-million-row
    file. If that becomes a real bottleneck, revisit cyvcf2 under WSL2
    rather than fighting the Windows build.

    This already streams line-by-line via gzip.open and never accumulates
    rows, so memory stays flat regardless of file size (a 2GB VCF is fine).
    The thing to watch is how you CONSUME it: `next(m for m in
    read_clinvar_vcf(path) if ...)` stays lazy and stops at the first
    match; `list(read_clinvar_vcf(path))` would force every one of
    ClinVar's 1M+ records into memory at once -- avoid that regardless of
    how much RAM you have.
    """
    import gzip

    opener = gzip.open if vcf_path.endswith(".gz") else open

    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue  # header/meta lines

            fields = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt, _qual, _filter, info_str = fields[:8]

            info = {}
            for item in info_str.split(";"):
                if "=" in item:
                    key, _, value = item.partition("=")
                    info[key] = value
                else:
                    info[item] = True  # flag field, no value

            gene_info = info.get("GENEINFO", "")          # e.g. "HBB:3043"
            gene = gene_info.split(":")[0] if gene_info else None

            # CLNHGVS often carries the transcript-level HGVS string, e.g.
            # "NM_000518.5(HBB):c.20A>T". Pull the transcript id out of it.
            # (Adjust if your ClinVar release names this field differently.)
            hgvs = info.get("CLNHGVS", "")
            transcript_id = hgvs.split("(")[0] if hgvs else ""

            # ClinVar rows are normally single-allele, but ALT can in
            # general be comma-separated -- split defensively.
            alt_allele = alt.split(",")[0]

            yield MutationRecord(
                gene=gene,
                transcript_id=transcript_id,
                chrom=chrom,
                genomic_position=int(pos),
                ref_allele=ref,
                alt_allele=alt_allele,
                clinical_significance=info.get("CLNSIG"),
                condition=info.get("CLNDN"),
            )


def read_clinvar_variant_summary(tsv_path: str, gene_symbol: Optional[str] = None, assembly: Optional[str] = "GRCh38"):
    """
    Stream-parse ClinVar's tab-delimited variant_summary.txt(.gz) one row
    at a time, yielding matching MutationRecord objects.

    NOT using pandas.read_csv() here anymore. variant_summary.txt has
    millions of rows across many string columns; pd.read_csv() reads the
    WHOLE file into a DataFrame before any filtering happens, and a
    DataFrame of that shape costs meaningfully more RAM than the raw file
    size (per-cell object overhead adds up fast at this scale) -- easily
    enough to exhaust memory on a multi-GB download, regardless of how
    narrow the eventual gene_symbol filter is. csv.DictReader + gzip
    processes one row at a time and only ever holds the current row (plus
    whatever MutationRecords you keep from it) in memory, independent of
    file size -- matching how read_clinvar_vcf already works above.

    Filter to one gene up front (gene_symbol="HBB") -- the full file covers
    every gene ClinVar has ever received a submission for, so without this
    you'd still be iterating everything even though memory stays flat.

    IMPORTANT: variant_summary.txt reports EVERY variant against BOTH
    GRCh37 and GRCh38 (one row each, distinguished by the Assembly column).
    Without filtering on assembly, you can pull a GRCh37 coordinate for a
    variant and use it against a GRCh38-parsed sequence -- wrong position,
    no error. Defaults to "GRCh38" to match ProtoDNASequencer's FASTA;
    pass assembly=None to disable the filter (get both).
    """
    import csv
    import gzip

    opener = gzip.open if tsv_path.endswith(".gz") else open

    with opener(tsv_path, "rt", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if assembly and row.get("Assembly") != assembly:
                continue
            if gene_symbol and gene_symbol not in (row.get("GeneSymbol") or "").split(";"):
                continue

            name = row.get("Name") or ""
            transcript_id = name.split("(")[0] if name else ""
            protein_change = None
            if "(p." in name:
                protein_change = "p." + name.split("(p.")[1].rstrip(")")

            # Not every row has a usable VCF-style position (some variant
            # types in this file -- large structural variants especially --
            # leave PositionVCF blank). Skip those rather than crashing the
            # whole stream on one malformed row.
            try:
                genomic_position = int(row["PositionVCF"])
            except (KeyError, ValueError):
                continue

            yield MutationRecord(
                gene=row.get("GeneSymbol"),
                transcript_id=transcript_id,
                chrom=row.get("Chromosome"),
                genomic_position=genomic_position,
                ref_allele=row.get("ReferenceAlleleVCF"),
                alt_allele=row.get("AlternateAlleleVCF"),
                protein_change=protein_change,
                clinical_significance=row.get("ClinicalSignificance"),
                condition=row.get("PhenotypeList"),
                assembly=row.get("Assembly"),
                review_status=row.get("ReviewStatus"),
            )


# ---------------------------------------------------------------------------
# 2. GTF/annotation side: modular lookup interface
# ---------------------------------------------------------------------------

class Strand(Enum):
    PLUS = "+"
    MINUS = "-"


@dataclass
class CodonLocation:
    """
    Genomic coordinates for ONE codon (3 bases) of a transcript's CDS.

    genomic_positions is always ascending (low to high) regardless of
    strand -- CodonComparer is responsible for reversing/complementing for
    minus-strand genes, not this class.
    """
    chrom: Optional[str]
    genomic_positions: tuple  # (int, int, int), ascending
    strand: Strand
    codon_number: int


class AnnotationLookup(ABC):
    """
    Abstract interface to an annotated-genome index (built from a GTF).
    PostgresAnnotationLookup (real) and MockAnnotationLookup (offline
    testing) both just need to satisfy this contract -- CodonComparer never
    needs to know which one it's talking to.
    """

    @abstractmethod
    def get_codon_location(self, transcript_id: str, codon_number: int) -> CodonLocation:
        """Genomic location of a given codon number in a transcript's CDS."""
        raise NotImplementedError

    @abstractmethod
    def get_codon_number_for_position(self, transcript_id: str, genomic_position: int) -> Optional[int]:
        """
        Which codon number (1-based) a genomic position falls into for a
        transcript, or None if it's outside that transcript's CDS
        (intronic/UTR).
        """
        raise NotImplementedError


class PostgresAnnotationLookup(AnnotationLookup):
    """
    Real AnnotationLookup, backed by the Postgres 'annotations' table via
    find_cds() (not find_exons() -- see below).

    IMPORTANT CAVEATS -- read before trusting this beyond a quick test:

    1. [RESOLVED] This now queries feature='CDS' rows via find_cds(),
       not feature='exon' rows. Exon boundaries include any 5'/3' UTR;
       CDS boundaries mark only the actual coding portion. Using exons
       here previously caused every codon number to be off by a constant
       amount equal to the UTR's length in codons -- confirmed concretely
       against HBB, which has a 50nt 5' UTR (NM_000518 CDS starts at
       mRNA position 51): codon 7 was coming back as codon 24, exactly
       17 codons (51nt) too high, matching the UTR length.

    2. A codon whose 3 bases straddle a splice junction (span two CDS
       blocks) is NOT handled -- get_codon_location raises
       NotImplementedError for that case rather than returning wrong
       coordinates. Extending CodonComparer._extract_codon to pull from
       two genomic ranges instead of one contiguous slice is required
       before this covers that case.

    3. 'chromosome' isn't available at the CDS/transcript level in the
       given schema (only find_gene/find_gene_by_name select seqid), so
       CodonLocation.chrom comes back None here. CodonComparer never reads
       it, so this isn't a functional blocker -- just flagging it's not
       populated.
    """

    def __init__(self, find_transcripts_fn=find_transcripts, find_cds_fn=find_cds):
        # Passed in as parameters (with these as defaults) rather than only
        # ever using the module-level import directly, so this class can be
        # unit-tested with fake functions instead of a live DB connection.
        self._find_transcripts = find_transcripts_fn
        self._find_cds = find_cds_fn
        self._cds_cache = {}   # transcript_id -> (ordered_cds_blocks, strand)

    def _ordered_cds(self, transcript_id):
        """
        CDS blocks for a transcript, in transcript (5'->3') order, each
        annotated with its cumulative CDS-relative offset. Memoized per
        transcript_id -- boundaries don't change between calls, and this
        avoids re-querying Postgres for every single codon lookup.
        """
        if transcript_id in self._cds_cache:
            return self._cds_cache[transcript_id]

        cds_blocks = self._find_cds(transcript_id)
        if not cds_blocks:
            raise KeyError(f"No CDS blocks found for transcript '{transcript_id}'.")

        strand = Strand(cds_blocks[0]["strand"])
        # find_cds() already returns rows ordered ascending by genomic
        # start. For a minus-strand transcript, 5'->3' runs in descending
        # genomic order, so reverse to get transcript order.
        ordered = list(cds_blocks) if strand is Strand.PLUS else list(reversed(cds_blocks))

        cumulative = 0
        for block in ordered:
            block_length = block["end"] - block["start"] + 1
            block["_cds_offset_start"] = cumulative   # 0-based, CDS-relative
            block["_length"] = block_length
            cumulative += block_length

        self._cds_cache[transcript_id] = (ordered, strand)
        return ordered, strand

    def get_codon_number_for_position(self, transcript_id: str, genomic_position: int) -> Optional[int]:
        ordered, strand = self._ordered_cds(transcript_id)

        for block in ordered:
            if block["start"] <= genomic_position <= block["end"]:
                if strand is Strand.PLUS:
                    within_block_offset = genomic_position - block["start"]
                else:
                    within_block_offset = block["end"] - genomic_position
                cds_offset = block["_cds_offset_start"] + within_block_offset
                return cds_offset // 3 + 1

        return None  # not in any CDS block of this transcript (intronic/UTR)

    def get_codon_location(self, transcript_id: str, codon_number: int) -> CodonLocation:
        ordered, strand = self._ordered_cds(transcript_id)

        cds_start = (codon_number - 1) * 3
        genomic_positions = [
            self._cds_offset_to_genomic(ordered, strand, cds_start + i)
            for i in range(3)
        ]

        if max(genomic_positions) - min(genomic_positions) != 2:
            raise NotImplementedError(
                f"Codon {codon_number} of {transcript_id} spans a splice "
                f"junction -- its 3 bases aren't 3 contiguous genomic "
                f"positions. CodonComparer/_extract_codon needs extending "
                f"to handle this case."
            )

        return CodonLocation(
            chrom=None,  # see class docstring, caveat 3
            genomic_positions=tuple(sorted(genomic_positions)),
            strand=strand,
            codon_number=codon_number,
        )

    @staticmethod
    def _cds_offset_to_genomic(ordered_cds, strand, cds_offset):
        for block in ordered_cds:
            if block["_cds_offset_start"] <= cds_offset < block["_cds_offset_start"] + block["_length"]:
                within_block_offset = cds_offset - block["_cds_offset_start"]
                if strand is Strand.PLUS:
                    return block["start"] + within_block_offset
                else:
                    return block["end"] - within_block_offset
        raise ValueError(f"CDS offset {cds_offset} is past the end of this transcript's CDS.")

class MockAnnotationLookup(AnnotationLookup):
    """
    Minimal in-memory stand-in for offline dev/tests that shouldn't need a
    live Postgres connection. Uses the REAL genomic coordinates for HBB
    codon 7 (GRCh38 chr11), derived as follows: HGVS c.20 is confirmed to
    be the exact middle nucleotide of codon 7, and reversing a 3-element
    sequence for a minus-strand gene always leaves the middle element in
    the middle -- so the codon spans 5,227,001-5,227,003.
    """

    def __init__(self):
        self._codon = CodonLocation(
            chrom="NC_000011.10",
            genomic_positions=(5_227_001, 5_227_002, 5_227_003),
            strand=Strand.MINUS,
            codon_number=7,
        )

    def get_codon_location(self, transcript_id: str, codon_number: int) -> CodonLocation:
        if transcript_id.startswith("NM_000518") and codon_number == 7:
            return self._codon
        raise KeyError(f"No mock data for {transcript_id} codon {codon_number}")

    def get_codon_number_for_position(self, transcript_id: str, genomic_position: int) -> Optional[int]:
        if transcript_id.startswith("NM_000518") and genomic_position in self._codon.genomic_positions:
            return self._codon.codon_number
        return None


# ---------------------------------------------------------------------------
# 3. The comparison itself
# ---------------------------------------------------------------------------

class CodonStatus(Enum):
    HEALTHY = "healthy"        # matches the reference codon
    MUTATED = "mutated"        # matches the known mutant codon
    UNEXPECTED = "unexpected"  # matches neither -- flag it, don't assume healthy


@dataclass
class CodonComparisonResult:
    gene: str
    codon_number: int
    reference_codon: str
    observed_codon: str
    mutant_codon: str
    reference_amino_acid: str
    observed_amino_acid: str
    status: CodonStatus


class CodonComparer:
    """
    Combines a parsed ChromosomeData, an AnnotationLookup, and a
    MutationRecord to answer: does the actual sequence at this codon match
    the healthy reference, the known mutant, or neither?

    NOTE: assumes a single-nucleotide substitution located within the
    codon (true for the sickle-cell example). ClinVar also contains indels
    and multi-base changes -- those need extra handling (allele length
    checks, possible frame shift) before this logic applies as-is.
    """

    def __init__(self, annotation_lookup: AnnotationLookup):
        self.annotation_lookup = annotation_lookup

    def _extract_codon(self, chrom_data, location: CodonLocation) -> str:
        """
        Pull 3 bases for a codon out of a parsed ChromosomeData, respecting
        strand.

        Plus strand: read directly off chrom_data.seq in ascending genomic
        order -- transcript order matches genomic order.

        Minus strand: transcript 5'->3' runs opposite to genomic order, so
        read chrom_data.complement (already computed once at sequencing
        time -- this is exactly why keeping it precomputed was worth it)
        across the same span, then reverse it into transcript order.
        """
        start, _, end = location.genomic_positions  # already ascending

        if location.strand is Strand.PLUS:
            return str(chrom_data.seq[start - 1:end]).upper()
        else:
            comp_slice_ascending = str(chrom_data.complement[start - 1:end]).upper()
            return comp_slice_ascending[::-1]

    def compare(self, chrom_data, mutation: MutationRecord) -> CodonComparisonResult:
        codon_number = self.annotation_lookup.get_codon_number_for_position(
            mutation.transcript_id, mutation.genomic_position
        )
        if codon_number is None:
            raise ValueError(
                f"Position {mutation.genomic_position} on "
                f"{mutation.transcript_id} isn't in a known CDS codon."
            )

        location = self.annotation_lookup.get_codon_location(mutation.transcript_id, codon_number)
        observed_codon = self._extract_codon(chrom_data, location)

        offset = self._offset_within_codon(location, mutation.genomic_position)
        reference_codon = self._substitute(observed_codon, offset, mutation.ref_allele, location.strand)
        mutant_codon = self._substitute(observed_codon, offset, mutation.alt_allele, location.strand)

        if observed_codon == reference_codon:
            status = CodonStatus.HEALTHY
        elif observed_codon == mutant_codon:
            status = CodonStatus.MUTATED
        else:
            status = CodonStatus.UNEXPECTED

        return CodonComparisonResult(
            gene=mutation.gene,
            codon_number=codon_number,
            reference_codon=reference_codon,
            observed_codon=observed_codon,
            mutant_codon=mutant_codon,
            reference_amino_acid=str(Seq(reference_codon).translate()),
            observed_amino_acid=str(Seq(observed_codon).translate()),
            status=status,
        )

    @staticmethod
    def _offset_within_codon(location: CodonLocation, genomic_position: int) -> int:
        """0/1/2 offset of a genomic position within the codon, in
        transcript (5'->3') order."""
        sorted_positions = sorted(location.genomic_positions)
        genome_index = sorted_positions.index(genomic_position)
        if location.strand is Strand.PLUS:
            return genome_index
        return 2 - genome_index  # minus strand: transcript order is reversed

    @staticmethod
    def _substitute(codon: str, offset: int, allele: str, strand: Strand) -> str:
        """Insert a ClinVar allele (always genome + strand) into a codon
        that's already in transcript orientation, complementing first if
        the gene is on the minus strand."""
        base = allele.upper()
        if strand is Strand.MINUS:
            base = str(Seq(base).complement())
        return codon[:offset] + base + codon[offset + 1:]


# ---------------------------------------------------------------------------
# Example usage
# ---------------------------------------------------------------------------

def main():
    p = ProtoDNASequencer("data/GCF_000001405.40_GRCh38.p14_genomic.fna")
    p.load_sequence()

    # sequenceChromosome(), not getChromosome() -- the sequencer is lazy now,
    # so chr11 has to be explicitly parsed before anything can read from it.
    chrom_data = p.sequenceChromosome("NC_000011.10")

    # Real DB-backed lookup. Swap for MockAnnotationLookup() if Postgres
    # isn't running / you just want to sanity-check the plumbing offline.
    lookup = PostgresAnnotationLookup()
    comparer = CodonComparer(lookup)

    # Streams the TSV row-by-row (see read_clinvar_variant_summary's
    # docstring) rather than loading the whole multi-GB file into memory.
    # gene_symbol="HBB" alone still returns every HBB variant ClinVar has
    # -- the position/allele filter below picks out the specific
    # sickle-cell record.
    mutations = read_clinvar_variant_summary(str(CLINVAR_VARIANT_SUMMARY), gene_symbol="HBB")
    mutation = next(
        m for m in mutations
        if m.genomic_position == 5_227_002 and m.ref_allele == "T" and m.alt_allele == "A"
    )

    result = comparer.compare(chrom_data, mutation)
    print(result)


if __name__ == "__main__":
    main()