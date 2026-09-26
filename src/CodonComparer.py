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
                             and strand. Your teammate's real, Postgres-
                             backed GTF index implements this; a minimal
                             in-memory stub is included so this file runs
                             standalone before that database exists.
3. CodonComparer          -> pulls the 3 bases at those coordinates out of
                             a parsed ChromosomeData (from the sequencer
                             project), applies ClinVar's ref->alt
                             substitution, and reports HEALTHY / MUTATED /
                             UNEXPECTED.

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

Install: pip install pandas biopython
(No VCF-parsing library required -- see read_clinvar_vcf's docstring for why.)
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum
from typing import Optional

from ProtoDNASequencer import ProtoDNASequencer

from Bio.Seq import Seq


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


def read_clinvar_variant_summary(tsv_path: str, gene_symbol: Optional[str] = None):
    """
    Parse ClinVar's tab-delimited variant_summary.txt(.gz) and yield
    MutationRecord objects. Includes HGVS protein-change strings the raw
    VCF doesn't carry on its own.

    Filter to one gene up front (gene_symbol="HBB") -- the full file covers
    every gene ClinVar has ever received a submission for.
    """
    import pandas as pd

    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    if gene_symbol:
        df = df[df["GeneSymbol"] == gene_symbol]

    for _, row in df.iterrows():
        name = row.get("Name", "")
        transcript_id = name.split("(")[0] if isinstance(name, str) else ""
        protein_change = None
        if isinstance(name, str) and "(p." in name:
            protein_change = "p." + name.split("(p.")[1].rstrip(")")

        yield MutationRecord(
            gene=row["GeneSymbol"],
            transcript_id=transcript_id,
            chrom=str(row["Chromosome"]),
            genomic_position=int(row["PositionVCF"]),
            ref_allele=row["ReferenceAlleleVCF"],
            alt_allele=row["AlternateAlleleVCF"],
            protein_change=protein_change,
            clinical_significance=row.get("ClinicalSignificance"),
            condition=row.get("PhenotypeList"),
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
    chrom: str
    genomic_positions: tuple  # (int, int, int), ascending
    strand: Strand
    codon_number: int


class AnnotationLookup(ABC):
    """
    Abstract interface to an annotated-genome index (built from a GTF).

    Your teammate's real, Postgres-backed implementation and any local
    test implementation (e.g. gffutils-based, see the skeleton at the
    bottom of this file) both just need to satisfy this contract.
    CodonComparer never needs to know which one it's talking to.
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


class MockAnnotationLookup(AnnotationLookup):
    """
    Minimal in-memory stand-in for local development, before the real
    Postgres-backed GTF index exists. Uses small, clearly-fake coordinates
    (NOT real GRCh38 positions) so the plumbing can be exercised without
    asserting unverified genomic coordinates as fact. Swap for the real
    implementation later -- nothing else in this file changes.
    """

    def __init__(self):
        self._codon = CodonLocation(
            chrom="11",
            genomic_positions=(100, 101, 102),   # toy positions, not real
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
        read chrom_data.complement (already computed once at load time --
        this is exactly why keeping it precomputed was worth it) across the
        same span, then reverse it into transcript order.
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
    # Swap in real data: chrom_data = your_sequencer.getChromosome("NC_000011.10")
    # This FakeChromData uses small toy coordinates matching MockAnnotationLookup's
    # toy positions (100-102) -- NOT real GRCh38 sequence.
    # class FakeChromData:
    #     seq = Seq("N" * 100 + "ACT" + "N" * 50)          # positions 101-103 -> "ACT"
    #     complement = Seq("N" * 100 + "TGA" + "N" * 50)    # complement of the above

    p = ProtoDNASequencer("data/GCF_000001405.40_GRCh38.p14_genomic.fna")
    p.load_sequence()
    chrom_data = p.getChromosome("NC_000011.10")

    lookup = MockAnnotationLookup()
    comparer = CodonComparer(lookup)

    mutation = MutationRecord(
        gene="HBB",
        transcript_id="NM_000518.5",
        chrom="11",
        genomic_position=101,
        ref_allele="A",
        alt_allele="C",
        protein_change="p.Glu7Val",
        clinical_significance="Pathogenic",
    )

    result = comparer.compare(chrom_data, mutation)
    print(result)


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------------
# Sketch: what a real AnnotationLookup implementation has to solve
# ---------------------------------------------------------------------------
# class GffutilsAnnotationLookup(AnnotationLookup):
#     """
#     Skeleton only -- the hard part isn't the library call, it's mapping a
#     codon number to genomic coordinates correctly when a transcript's CDS
#     spans multiple exons (a codon can straddle a splice junction, so its
#     3 bases aren't always 3 consecutive genomic positions).
#
#     General approach:
#       1. Get the transcript's CDS features, ordered 5'->3' along the
#          transcript (reverse genomic order for minus-strand genes).
#       2. Build a cumulative-length index across those CDS chunks.
#       3. codon_number -> cDNA offset (codon_number - 1) * 3 -> walk the
#          cumulative index to find which CDS chunk(s) that offset falls
#          in, and translate back to genomic coordinate(s) -- possibly 2
#          genomic ranges if the codon straddles a splice junction.
#     Your teammate's Postgres schema may already store CDS-relative
#     coordinates precomputed, which would make this much simpler than
#     doing it here -- worth checking before reimplementing it.
#     """
#     def __init__(self, gtf_path: str, db_path: str = ":memory:"):
#         import gffutils
#         self.db = gffutils.create_db(gtf_path, db_path, force=True, keep_order=True)
#
#     def get_codon_location(self, transcript_id, codon_number):
#         raise NotImplementedError("see docstring above")
#
#     def get_codon_number_for_position(self, transcript_id, genomic_position):
#         raise NotImplementedError("see docstring above")