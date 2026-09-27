import sys
from pathlib import Path

import psycopg
from psycopg.rows import dict_row

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(PROJECT_ROOT))

from src.ProtoDNASequencer import ProtoDNASequencer


FASTA_FILE = PROJECT_ROOT / "GCF_000001405.40_GRCh38.p14_genomic.fna"

# sequencer = ProtoDNASequencer(FASTA_FILE)
# sequencer.load_sequence()


# ============================================================
# DATABASE CONFIG
# Must match the config used by your importer script.
# ============================================================

DB_CONFIG = {
    "host": "localhost",
    "port": 5432,
    "dbname": "genome",
    "user": "postgres",
    "password": "password",
}


def get_connection():
    """
    Opens a new connection to the annotations database.

    Rows are returned as dictionaries (dict_row) so the rest of this
    file can keep returning plain dicts, just like the old GTF-based
    version did.
    """
    return psycopg.connect(**DB_CONFIG, row_factory=dict_row)


# ============================================================
# LOOKUPS
# ============================================================

def find_gene_by_name(gene_name):

    query = """
        SELECT
            seqid AS chromosome,
            gene_name AS gene,
            gene_id,
            start_position AS start,
            end_position AS "end",
            strand
        FROM annotations
        WHERE feature = 'gene'
          AND (gene_name = %(gene_name)s
               OR attributes ->> 'gene_name' = %(gene_name)s)
        LIMIT 1;
    """

    with get_connection() as connection:
        with connection.cursor() as cursor:
            cursor.execute(query, {"gene_name": gene_name})
            return cursor.fetchone()


def find_gene(chromosome, position):

    query = """
        SELECT
            gene_name AS gene,
            seqid AS chromosome,
            gene_id,
            start_position AS start,
            end_position AS "end",
            strand
        FROM annotations
        WHERE feature = 'gene'
          AND seqid = %(chromosome)s
          AND start_position <= %(position)s
          AND end_position >= %(position)s
        LIMIT 1;
    """

    with get_connection() as connection:
        with connection.cursor() as cursor:
            cursor.execute(
                query,
                {"chromosome": chromosome, "position": position},
            )
            return cursor.fetchone()


def find_transcripts(gene_id):

    query = """
        SELECT
            transcript_id,
            attributes ->> 'transcript_name' AS transcript_name,
            gene_id,
            start_position AS start,
            end_position AS "end",
            strand
        FROM annotations
        WHERE feature = 'transcript'
          AND gene_id = %(gene_id)s
        ORDER BY start_position;
    """

    with get_connection() as connection:
        with connection.cursor() as cursor:
            cursor.execute(query, {"gene_id": gene_id})
            return cursor.fetchall()


def find_exons(transcript_id):

    query = """
        SELECT
            attributes ->> 'exon_id' AS exon_id,
            transcript_id,
            start_position AS start,
            end_position AS "end",
            strand
        FROM annotations
        WHERE feature = 'exon'
          AND transcript_id = %(transcript_id)s
        ORDER BY start_position;
    """

    with get_connection() as connection:
        with connection.cursor() as cursor:
            cursor.execute(query, {"transcript_id": transcript_id})
            return cursor.fetchall()


def find_cds(transcript_id):
    """
    CDS (coding sequence) blocks for a transcript -- NOT the same as
    find_exons(). Exon boundaries include any 5'/3' UTR; CDS boundaries
    mark only the actual coding portion, which is what codon-number math
    needs. For a transcript with a 5' UTR, using find_exons() here would
    systematically offset every codon number by the UTR's length (this is
    exactly the bug that showed up testing HBB codon 7 -- it does have a
    50nt 5' UTR, and every codon number came back 17 too high as a result).
    """

    query = """
        SELECT
            transcript_id,
            start_position AS start,
            end_position AS "end",
            strand
        FROM annotations
        WHERE feature = 'CDS'
          AND transcript_id = %(transcript_id)s
        ORDER BY start_position;
    """

    with get_connection() as connection:
        with connection.cursor() as cursor:
            cursor.execute(query, {"transcript_id": transcript_id})
            return cursor.fetchall()


# ============================================================
# SEQUENCE LOOKUP (unchanged - still reads from the FASTA file
# via ProtoDNASequencer, since the database only stores
# annotations, not the raw sequence itself)
# ============================================================


# ============================================================
# EXAMPLE USAGE
# ============================================================

if __name__ == "__main__":

    first_chromosome_gene = find_gene_by_name("HBB")
    print(first_chromosome_gene)

    if first_chromosome_gene:
        print(
            first_chromosome_gene["chromosome"],
            first_chromosome_gene["start"],
            first_chromosome_gene["end"],
        )

        transcripts = find_transcripts(first_chromosome_gene["gene_id"])
        print(transcripts)

        if transcripts:
            exons = find_exons(transcripts[0]["transcript_id"])
            print(exons)

            if exons:
                print(
                    first_chromosome_gene["chromosome"],
                    exons[0]["start"],
                    exons[0]["end"],
                )

            cds = find_cds(transcripts[0]["transcript_id"])
            print(cds)