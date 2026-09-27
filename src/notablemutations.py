"""
Scan ClinVar's variant_summary.txt(.gz) ONCE for a curated list of
"identifiable" single-nucleotide, high-confidence pathogenic mutations,
and cache the result as a small JSON file. Your UI loads that small file
instantly instead of re-scanning the full multi-GB file every run.

Filters applied (each one exists to cut a specific kind of noise out of
a raw ClinVar dump):

- Assembly == "GRCh38" -- variant_summary.txt reports every variant
  against BOTH GRCh37 and GRCh38. Skipping this filter mixes coordinate
  systems; see read_clinvar_variant_summary's docstring in CodonComparer.py.

- ClinicalSignificance STARTS WITH "pathogenic" or "likely pathogenic"
  -- not a substring check. "Conflicting classifications of pathogenicity"
  contains the substring "pathogenic" too, but means ClinVar submitters
  disagree, which is the opposite of what "notable, well-vetted" means
  here.

- ReviewStatus is one of ClinVar's two highest-confidence tiers:
  "reviewed by expert panel" or "practice guideline". This is what
  actually narrows millions of rows down to a genuinely curated subset --
  most ClinVar rows are single-submitter assertions with no independent
  review.

- Single-nucleotide substitution only (len(ref) == 1 and len(alt) == 1)
  -- this is what "identifiable single [codon] change" means concretely.
  Excludes the large deletions/duplications/frameshift indels that
  dominate a raw ClinVar dump for any given gene (see the HBB dump from
  earlier -- most of it was exactly this kind of noise).

Only chromosomes in GRCH38_CHROM_ACCESSIONS (the 22 autosomes + X/Y/MT)
are kept -- anything on an unplaced scaffold is skipped, since your
ProtoDNASequencer's chromosomes_only=True filtering wouldn't have that
sequence loaded anyway.

Run once -- it's still a full pass over the file, so expect this to take
a few minutes even though memory stays flat (see read_clinvar_variant_summary's
streaming implementation). Re-run only when you want to refresh the cache
against a newer ClinVar release.
"""

import json
from pathlib import Path

from CodonComparer import (
    read_clinvar_variant_summary,
    CLINVAR_VARIANT_SUMMARY,
    GRCH38_CHROM_ACCESSIONS,
)

HIGH_CONFIDENCE_REVIEW_STATUSES = {"reviewed by expert panel", "practice guideline"}


def is_pathogenic(significance):
    """True for a confirmed pathogenic/likely-pathogenic call -- deliberately
    NOT a substring check (see module docstring: "Conflicting classifications
    of pathogenicity" would otherwise false-match)."""
    if not significance:
        return False
    sig = significance.lower()
    return sig.startswith("pathogenic") or sig.startswith("likely pathogenic")


def find_notable_mutations(tsv_path, output_path="notable_mutations.json", progress_every=500_000):
    """
    Streams tsv_path once, applying the filters described above, and
    writes matches to output_path as a JSON list. Returns that list too.
    """
    results = []
    seen = set()  # (gene, position, ref, alt) -- guards against exact duplicate rows

    mutations = read_clinvar_variant_summary(tsv_path, gene_symbol=None, assembly="GRCh38")

    scanned = 0
    for m in mutations:
        scanned += 1
        if progress_every and scanned % progress_every == 0:
            print(f"...scanned {scanned:,} rows, {len(results)} matches so far")

        if not m.ref_allele or not m.alt_allele:
            continue
        if len(m.ref_allele) != 1 or len(m.alt_allele) != 1:
            continue  # not a single-nucleotide substitution

        if not is_pathogenic(m.clinical_significance):
            continue
        if (m.review_status or "").lower() not in HIGH_CONFIDENCE_REVIEW_STATUSES:
            continue

        chrom_accession = GRCH38_CHROM_ACCESSIONS.get(m.chrom)
        if chrom_accession is None:
            continue  # unplaced scaffold or unrecognized chromosome value

        key = (m.gene, m.genomic_position, m.ref_allele, m.alt_allele)
        if key in seen:
            continue
        seen.add(key)

        results.append({
            "gene": m.gene,
            "transcript_id": m.transcript_id,
            "chrom_accession": chrom_accession,
            "genomic_position": m.genomic_position,
            "ref_allele": m.ref_allele,
            "alt_allele": m.alt_allele,
            "protein_change": m.protein_change,
            "clinical_significance": m.clinical_significance,
            "review_status": m.review_status,
            "condition": m.condition,
        })

    Path(output_path).write_text(json.dumps(results, indent=2))
    print(f"Scanned {scanned:,} rows total. Wrote {len(results)} notable mutations to {output_path}")
    return results


if __name__ == "__main__":
    find_notable_mutations(str(CLINVAR_VARIANT_SUMMARY))