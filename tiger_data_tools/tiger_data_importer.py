import psycopg
from psycopg import sql

from pathlib import Path
import json
import re
from pathlib import Path

import psycopg
from psycopg.types.json import Jsonb


# ============================================================
# CONFIGURATION
# ============================================================

DB_CONFIG = {
    "host": "localhost",
    "port": 5432,
    "dbname": "genome",
    "user": "postgres",
    "password": "password",
}

# CHANGE THIS to the location of your GTF file.
GTF_FILE = Path(r"data/genomic.gtf")

# How many rows to send in one batch.
BATCH_SIZE = 10_000


# ============================================================
# ATTRIBUTE PARSING
# ============================================================

def add_attribute(attributes, key, value):
    """
    Adds an attribute while preserving duplicate keys.

    Example:
        db_xref "GeneID:123";
        db_xref "HGNC:456";

    becomes:

        {
            "db_xref": [
                "GeneID:123",
                "HGNC:456"
            ]
        }
    """
    if key in attributes:
        if isinstance(attributes[key], list):
            attributes[key].append(value)
        else:
            attributes[key] = [attributes[key], value]
    else:
        attributes[key] = value


def parse_attributes(raw_attributes):
    """
    Parse a GTF/GFF3 column-9 attribute string.

    Supports normal GTF:
        gene_id "123"; gene_name "HBB";

    and GFF3-style:
        ID=gene1;Parent=gene0;

    Duplicate attributes are preserved as lists.
    """

    attributes = {}

    # Remove trailing whitespace/semicolon.
    text = raw_attributes.strip().rstrip(";").strip()

    if not text:
        return attributes

    # --------------------------------------------------------
    # GTF style: key "value";
    # Handles escaped quotes inside quoted values.
    # --------------------------------------------------------
    gtf_pattern = re.compile(
        r'([^\s;=]+)\s+"((?:\\.|[^"\\])*)"\s*(?:;|$)'
    )

    gtf_matches = list(gtf_pattern.finditer(text))

    if gtf_matches:
        for match in gtf_matches:
            key = match.group(1)
            value = match.group(2)

            # Unescape common GTF escapes.
            value = value.replace(r'\"', '"').replace(r'\\', '\\')

            add_attribute(attributes, key, value)

        # There can occasionally be unquoted attributes mixed in.
        consumed = [False] * len(text)
        for match in gtf_matches:
            for i in range(match.start(), match.end()):
                consumed[i] = True

        leftover = "".join(
            char if not consumed[i] else " "
            for i, char in enumerate(text)
        )

        # Parse anything left over as GFF3/unquoted syntax.
        for piece in leftover.split(";"):
            piece = piece.strip()

            if not piece:
                continue

            if "=" in piece:
                key, value = piece.split("=", 1)
                add_attribute(
                    attributes,
                    key.strip(),
                    value.strip().strip('"')
                )

        return attributes

    # --------------------------------------------------------
    # GFF3/unquoted style: key=value;
    # --------------------------------------------------------
    for piece in text.split(";"):
        piece = piece.strip()

        if not piece:
            continue

        if "=" in piece:
            key, value = piece.split("=", 1)
            key = key.strip()
            value = value.strip()

            # GFF3 can use URL-style escaping.
            value = (
                value
                .replace("%3B", ";")
                .replace("%3D", "=")
                .replace("%2C", ",")
                .replace("%25", "%")
            )

            add_attribute(attributes, key, value)

        else:
            # Last-resort handling for an unusual attribute.
            parts = piece.split(None, 1)

            if len(parts) == 2:
                key, value = parts
                add_attribute(
                    attributes,
                    key.strip(),
                    value.strip().strip('"')
                )

    return attributes


# ============================================================
# DATABASE SETUP
# ============================================================

def create_database_tables(connection):
    """
    Creates a generic annotation table.

    We intentionally do NOT create separate gene/transcript/exon/CDS
    tables. The feature column tells us what each row represents, so
    this works for every feature type in the file.
    """

    with connection.cursor() as cursor:

        cursor.execute("""
            DROP TABLE IF EXISTS annotations;
        """)

        cursor.execute("""
            DROP TABLE IF EXISTS gtf_metadata;
        """)

        # ----------------------------------------------------
        # Metadata / header lines
        # ----------------------------------------------------
        cursor.execute("""
            CREATE TABLE gtf_metadata (
                metadata_id BIGSERIAL PRIMARY KEY,
                line_number BIGINT NOT NULL,
                raw_line TEXT NOT NULL
            );
        """)

        # ----------------------------------------------------
        # Main annotation table
        # ----------------------------------------------------
        cursor.execute("""
            CREATE TABLE annotations (
                annotation_id BIGSERIAL PRIMARY KEY,

                -- Original GTF/GFF columns
                seqid TEXT NOT NULL,
                source TEXT NOT NULL,
                feature TEXT NOT NULL,
                start_position BIGINT NOT NULL,
                end_position BIGINT NOT NULL,
                score TEXT NOT NULL,
                strand CHAR(1) NOT NULL,
                frame CHAR(1) NOT NULL,

                -- Column 9 exactly as it appeared in the file
                attributes_raw TEXT NOT NULL,

                -- Parsed version of EVERY attribute
                attributes JSONB NOT NULL,

                -- Frequently searched attributes extracted into columns
                gene_id TEXT,
                transcript_id TEXT,
                gene_name TEXT,
                exon_number TEXT,
                parent_id TEXT,
                feature_id TEXT,
                protein_id TEXT
            );
        """)

    connection.commit()


# ============================================================
# INSERT BATCH
# ============================================================

def insert_batch(connection, rows):
    """
    Insert one batch of annotation rows.

    executemany is used here because it is simple and reliable.
    PostgreSQL indexes are created after the import so the initial
    load is substantially faster.
    """

    if not rows:
        return

    sql = """
        INSERT INTO annotations (
            seqid,
            source,
            feature,
            start_position,
            end_position,
            score,
            strand,
            frame,
            attributes_raw,
            attributes,
            gene_id,
            transcript_id,
            gene_name,
            exon_number,
            parent_id,
            feature_id,
            protein_id
        )
        VALUES (
            %s, %s, %s, %s, %s, %s, %s, %s,
            %s, %s, %s, %s, %s, %s, %s, %s, %s
        )
    """

    with connection.cursor() as cursor:
        cursor.executemany(sql, rows)

    connection.commit()


# ============================================================
# INDEXES
# ============================================================

def create_indexes(connection):
    """
    Creates indexes after the data has been imported.
    """

    print("Creating indexes...")

    with connection.cursor() as cursor:

        cursor.execute("""
            CREATE INDEX idx_annotations_feature
            ON annotations(feature);
        """)

        cursor.execute("""
            CREATE INDEX idx_annotations_seqid
            ON annotations(seqid);
        """)

        cursor.execute("""
            CREATE INDEX idx_annotations_gene_id
            ON annotations(gene_id);
        """)

        cursor.execute("""
            CREATE INDEX idx_annotations_transcript_id
            ON annotations(transcript_id);
        """)

        cursor.execute("""
            CREATE INDEX idx_annotations_gene_name
            ON annotations(gene_name);
        """)

        cursor.execute("""
            CREATE INDEX idx_annotations_parent_id
            ON annotations(parent_id);
        """)

        cursor.execute("""
            CREATE INDEX idx_annotations_feature_id
            ON annotations(feature_id);
        """)

        cursor.execute("""
            CREATE INDEX idx_annotations_protein_id
            ON annotations(protein_id);
        """)

        # Useful for queries such as:
        # "Find all genes on this chromosome."
        cursor.execute("""
            CREATE INDEX idx_annotations_seqid_feature
            ON annotations(seqid, feature);
        """)

        # Useful for searching arbitrary JSON attributes.
        cursor.execute("""
            CREATE INDEX idx_annotations_attributes
            ON annotations USING GIN(attributes);
        """)

        # Useful for genomic-region queries.
        cursor.execute("""
            CREATE INDEX idx_annotations_coordinates
            ON annotations(seqid, start_position, end_position);
        """)

    connection.commit()


# ============================================================
# GTF IMPORT
# ============================================================

def import_gtf(connection, gtf_path):
    """
    Reads the entire GTF/GFF3 file and inserts every annotation.
    """

    print(f"Reading: {gtf_path}")

    if not gtf_path.exists():
        raise FileNotFoundError(
            f"GTF file was not found:\n{gtf_path}"
        )

    rows = []
    metadata_rows = []

    annotation_count = 0
    metadata_count = 0

    with gtf_path.open(
        "r",
        encoding="utf-8",
        errors="replace"
    ) as file:

        for line_number, raw_line in enumerate(file, start=1):

            line = raw_line.rstrip("\n\r")

            if not line:
                continue

            # ------------------------------------------------
            # Header/directive line
            # ------------------------------------------------
            if line.startswith("#"):
                metadata_rows.append(
                    (line_number, line)
                )
                metadata_count += 1
                continue

            # ------------------------------------------------
            # GTF should contain exactly 9 tab-separated fields
            # ------------------------------------------------
            fields = line.split("\t")

            if len(fields) != 9:
                print(
                    f"WARNING: line {line_number} has "
                    f"{len(fields)} columns instead of 9. "
                    f"Skipping."
                )
                continue

            (
                seqid,
                source,
                feature,
                start,
                end,
                score,
                strand,
                frame,
                raw_attributes,
            ) = fields

            try:
                start = int(start)
                end = int(end)
            except ValueError:
                print(
                    f"WARNING: invalid coordinates on "
                    f"line {line_number}. Skipping."
                )
                continue

            attributes = parse_attributes(raw_attributes)

            # ------------------------------------------------
            # Extract commonly used values.
            #
            # NCBI's GTF uses:
            #   gene
            #   gene_biotype
            #   transcript_id
            #   transcript_biotype
            #   exon_number
            #
            # The COMPLETE attributes dict is still retained.
            # ------------------------------------------------
            gene_id = attributes.get("gene_id")
            transcript_id = attributes.get("transcript_id")
            gene_name = attributes.get("gene")
            exon_number = attributes.get("exon_number")
            parent_id = attributes.get("Parent")
            feature_id = attributes.get("ID")
            protein_id = attributes.get("protein_id")

            # If a duplicate attribute became a list, the database
            # column receives the JSON representation rather than
            # silently throwing information away.
            def scalar_or_json(value):
                if isinstance(value, list):
                    return json.dumps(value)
                return value

            rows.append((
                seqid,
                source,
                feature,
                start,
                end,
                score,
                strand,
                frame,
                raw_attributes,
                Jsonb(attributes),
                scalar_or_json(gene_id),
                scalar_or_json(transcript_id),
                scalar_or_json(gene_name),
                scalar_or_json(exon_number),
                scalar_or_json(parent_id),
                scalar_or_json(feature_id),
                scalar_or_json(protein_id),
            ))

            annotation_count += 1

            if len(rows) >= BATCH_SIZE:
                insert_batch(connection, rows)
                print(
                    f"Imported {annotation_count:,} "
                    f"annotation rows..."
                )
                rows.clear()

    # Insert remaining rows.
    if rows:
        insert_batch(connection, rows)

    # --------------------------------------------------------
    # Save header/directive information.
    # --------------------------------------------------------
    if metadata_rows:
        with connection.cursor() as cursor:
            cursor.executemany(
                """
                INSERT INTO gtf_metadata (
                    line_number,
                    raw_line
                )
                VALUES (%s, %s)
                """,
                metadata_rows
            )

        connection.commit()

    print()
    print("Import complete.")
    print(f"Annotation rows: {annotation_count:,}")
    print(f"Metadata lines:   {metadata_count:,}")


# ============================================================
# VERIFY IMPORT
# ============================================================

def verify_import(connection):
    """
    Runs a few sanity checks after import.
    """

    print()
    print("Verifying database...")

    with connection.cursor() as cursor:

        cursor.execute("""
            SELECT COUNT(*)
            FROM annotations;
        """)
        total = cursor.fetchone()[0]

        cursor.execute("""
            SELECT feature, COUNT(*)
            FROM annotations
            GROUP BY feature
            ORDER BY COUNT(*) DESC;
        """)
        feature_counts = cursor.fetchall()

        cursor.execute("""
            SELECT COUNT(*)
            FROM gtf_metadata;
        """)
        metadata_count = cursor.fetchone()[0]

    print(f"Total annotation rows: {total:,}")
    print(f"Metadata rows:          {metadata_count:,}")

    print()
    print("Feature types:")

    for feature, count in feature_counts:
        print(f"  {feature:<25} {count:,}")


# ============================================================
# MAIN
# ============================================================

def main():

    print("=" * 60)
    print("NCBI GTF -> PostgreSQL Importer")
    print("=" * 60)

    print()
    print(f"Database: {DB_CONFIG['dbname']}")
    print(f"GTF:      {GTF_FILE}")
    print()

    try:
        with psycopg.connect(**DB_CONFIG) as connection:

            print("Connected to PostgreSQL.")

            print("Creating tables...")
            create_database_tables(connection)

            import_gtf(
                connection,
                GTF_FILE
            )

            create_indexes(connection)

            verify_import(connection)

        print()
        print("=" * 60)
        print("DONE")
        print("=" * 60)

    except psycopg.OperationalError as error:
        print()
        print("Could not connect to PostgreSQL.")
        print()
        print(error)
        print()
        print("Check:")
        print("  1. PostgreSQL is running.")
        print("  2. Database 'genome' exists.")
        print("  3. Username/password are correct.")
        print("  4. PostgreSQL is using port 5432.")

    except Exception as error:
        print()
        print("IMPORT FAILED")
        print()
        print(type(error).__name__)
        print(error)
        raise


if __name__ == "__main__":
    main()


# path = Path("/mnt/data/postgres_importer.py")
# path.write_text(code, encoding="utf-8")

# print(f"Created: {path}")

# # Change these if you/used different PostgreSQL credentials/database.
# DB_CONFIG = {
#     "dbname": "genome",
#     "user": "postgres",
#     "password": "password",
#     "host": "localhost",
#     "port": 5432
# }


# def parse_attributes(attribute_string):
#     """
#     Convert the GTF attributes column into a dictionary.

#     Example:
#         gene_id "ENSG123"; gene_name "HBB";

#     becomes:
#         {
#             "gene_id": "ENSG123",
#             "gene_name": "HBB"
#         }
#     """

#     attributes = {}

#     for item in attribute_string.strip().split(";"):
#         item = item.strip()

#         if not item:
#             continue

#         key, value = item.split(" ", 1)
#         attributes[key] = value.strip('"')

#     return attributes


# def create_tables(conn):
#     """Create the tables used to store the GTF data."""

#     with conn.cursor() as cur:

#         cur.execute("""
#             CREATE TABLE IF NOT EXISTS genes (
#                 id BIGSERIAL PRIMARY KEY,
#                 gene_id TEXT,
#                 gene_name TEXT,
#                 gene_type TEXT,
#                 chromosome TEXT,
#                 start_position INTEGER,
#                 end_position INTEGER,
#                 score TEXT,
#                 strand CHAR(1),
#                 frame TEXT
#             );
#         """)

#         cur.execute("""
#             CREATE TABLE IF NOT EXISTS transcripts (
#                 id BIGSERIAL PRIMARY KEY,
#                 transcript_id TEXT,
#                 gene_id TEXT,
#                 gene_name TEXT,
#                 transcript_type TEXT,
#                 chromosome TEXT,
#                 start_position INTEGER,
#                 end_position INTEGER,
#                 score TEXT,
#                 strand CHAR(1),
#                 frame TEXT
#             );
#         """)

#         cur.execute("""
#             CREATE TABLE IF NOT EXISTS exons (
#                 id BIGSERIAL PRIMARY KEY,
#                 exon_id TEXT,
#                 transcript_id TEXT,
#                 gene_id TEXT,
#                 gene_name TEXT,
#                 chromosome TEXT,
#                 start_position INTEGER,
#                 end_position INTEGER,
#                 score TEXT,
#                 strand CHAR(1),
#                 frame TEXT
#             );
#         """)

#         cur.execute("""
#             CREATE TABLE IF NOT EXISTS cds (
#                 id BIGSERIAL PRIMARY KEY,
#                 protein_id TEXT,
#                 transcript_id TEXT,
#                 gene_id TEXT,
#                 gene_name TEXT,
#                 chromosome TEXT,
#                 start_position INTEGER,
#                 end_position INTEGER,
#                 score TEXT,
#                 strand CHAR(1),
#                 frame TEXT
#             );
#         """)

#         # Indexes for the queries your application will make frequently.
#         cur.execute("""
#             CREATE INDEX IF NOT EXISTS idx_genes_gene_name
#             ON genes(gene_name);
#         """)

#         cur.execute("""
#             CREATE INDEX IF NOT EXISTS idx_genes_gene_id
#             ON genes(gene_id);
#         """)

#         cur.execute("""
#             CREATE INDEX IF NOT EXISTS idx_transcripts_gene_id
#             ON transcripts(gene_id);
#         """)

#         cur.execute("""
#             CREATE INDEX IF NOT EXISTS idx_transcripts_transcript_id
#             ON transcripts(transcript_id);
#         """)

#         cur.execute("""
#             CREATE INDEX IF NOT EXISTS idx_exons_transcript_id
#             ON exons(transcript_id);
#         """)

#         cur.execute("""
#             CREATE INDEX IF NOT EXISTS idx_exons_gene_id
#             ON exons(gene_id);
#         """)

#         cur.execute("""
#             CREATE INDEX IF NOT EXISTS idx_cds_transcript_id
#             ON cds(transcript_id);
#         """)

#         cur.execute("""
#             CREATE INDEX IF NOT EXISTS idx_cds_gene_id
#             ON cds(gene_id);
#         """)

#     conn.commit()


# def import_gtf(conn, filename):

#     # Keep separate batches so we can efficiently insert many rows at once.
#     gene_batch = []
#     transcript_batch = []
#     exon_batch = []
#     cds_batch = []

#     BATCH_SIZE = 5000

#     with open(filename, "r", encoding="utf-8") as file:

#         for line_number, line in enumerate(file, start=1):

#             if line.startswith("#"):
#                 continue

#             fields = line.rstrip("\n").split("\t")

#             if len(fields) != 9:
#                 print(f"Skipping malformed line {line_number}")
#                 continue

#             chromosome = fields[0]
#             feature = fields[2]
#             start = int(fields[3])
#             end = int(fields[4])
#             score = fields[5]
#             strand = fields[6]
#             frame = fields[7]

#             attributes = parse_attributes(fields[8])

#             gene_id = attributes.get("gene_id")
#             gene_name = attributes.get("gene_name")

#             if feature == "gene":

#                 gene_batch.append((
#                     gene_id,
#                     gene_name,
#                     attributes.get("gene_type"),
#                     chromosome,
#                     start,
#                     end,
#                     score,
#                     strand,
#                     frame
#                 ))

#             elif feature == "transcript":

#                 transcript_batch.append((
#                     attributes.get("transcript_id"),
#                     gene_id,
#                     gene_name,
#                     attributes.get("transcript_type"),
#                     chromosome,
#                     start,
#                     end,
#                     score,
#                     strand,
#                     frame
#                 ))

#             elif feature == "exon":

#                 exon_batch.append((
#                     attributes.get("exon_id"),
#                     attributes.get("transcript_id"),
#                     gene_id,
#                     gene_name,
#                     chromosome,
#                     start,
#                     end,
#                     score,
#                     strand,
#                     frame
#                 ))

#             elif feature == "CDS":

#                 cds_batch.append((
#                     attributes.get("protein_id"),
#                     attributes.get("transcript_id"),
#                     gene_id,
#                     gene_name,
#                     chromosome,
#                     start,
#                     end,
#                     score,
#                     strand,
#                     frame
#                 ))

#             # Insert when a batch gets large enough.
#             if len(gene_batch) >= BATCH_SIZE:
#                 insert_genes(conn, gene_batch)
#                 gene_batch.clear()

#             if len(transcript_batch) >= BATCH_SIZE:
#                 insert_transcripts(conn, transcript_batch)
#                 transcript_batch.clear()

#             if len(exon_batch) >= BATCH_SIZE:
#                 insert_exons(conn, exon_batch)
#                 exon_batch.clear()

#             if len(cds_batch) >= BATCH_SIZE:
#                 insert_cds(conn, cds_batch)
#                 cds_batch.clear()

#             if line_number % 100000 == 0:
#                 print(f"Processed {line_number:,} lines...")

#     # Insert anything remaining.
#     if gene_batch:
#         insert_genes(conn, gene_batch)

#     if transcript_batch:
#         insert_transcripts(conn, transcript_batch)

#     if exon_batch:
#         insert_exons(conn, exon_batch)

#     if cds_batch:
#         insert_cds(conn, cds_batch)

#     print("GTF import complete.")


# def insert_genes(conn, rows):

#     with conn.cursor() as cur:
#         cur.executemany("""
#             INSERT INTO genes (
#                 gene_id,
#                 gene_name,
#                 gene_type,
#                 chromosome,
#                 start_position,
#                 end_position,
#                 score,
#                 strand,
#                 frame
#             )
#             VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
#         """, rows)

#     conn.commit()


# def insert_transcripts(conn, rows):

#     with conn.cursor() as cur:
#         cur.executemany("""
#             INSERT INTO transcripts (
#                 transcript_id,
#                 gene_id,
#                 gene_name,
#                 transcript_type,
#                 chromosome,
#                 start_position,
#                 end_position,
#                 score,
#                 strand,
#                 frame
#             )
#             VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
#         """, rows)

#     conn.commit()


# def insert_exons(conn, rows):

#     with conn.cursor() as cur:
#         cur.executemany("""
#             INSERT INTO exons (
#                 exon_id,
#                 transcript_id,
#                 gene_id,
#                 gene_name,
#                 chromosome,
#                 start_position,
#                 end_position,
#                 score,
#                 strand,
#                 frame
#             )
#             VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
#         """, rows)

#     conn.commit()


# def insert_cds(conn, rows):

#     with conn.cursor() as cur:
#         cur.executemany("""
#             INSERT INTO cds (
#                 protein_id,
#                 transcript_id,
#                 gene_id,
#                 gene_name,
#                 chromosome,
#                 start_position,
#                 end_position,
#                 score,
#                 strand,
#                 frame
#             )
#             VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
#         """, rows)

#     conn.commit()


# if __name__ == "__main__":

#     print("Connecting to PostgreSQL...")

#     with psycopg.connect(**DB_CONFIG) as conn:

#         print("Connected.")

#         create_tables(conn)

#         print("Importing GTF...")
#         import_gtf(conn, GTF_FILE)

#     print("Done.")
# # ```

# # ### Before running it

# # You need a PostgreSQL database called `genome`.

# # Once PostgreSQL is installed, you can create it with:

# # ```sql
# # CREATE DATABASE genome;
# # ```

# # Then change:

# # ```python
# # "password": "YOUR_PASSWORD"
# # ```

# # to the password you gave your PostgreSQL `postgres` user.

# # Your resulting structure will be:

# # ```text
# # genome
# # │
# # ├── genes
# # │    ├── HBB
# # │    ├── WASH7P
# # │    └── ...
# # │
# # ├── transcripts
# # │
# # ├── exons
# # │
# # └── cds
# # ```

# # And finding HBB becomes:

# # ```sql
# # SELECT *
# # FROM genes
# # WHERE gene_name = 'HBB';
# # ```

# # This is a better foundation for the coordinate-based queries you'll eventually need for the DNA visualizer.
