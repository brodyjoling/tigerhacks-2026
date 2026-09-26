# PostgreSQL Genome Database Setup

This guide explains how to set up the PostgreSQL genome annotation database locally and import the NCBI GTF file.

The database stores the genome annotation data used by our DNA visualization project.

---

## What This Does

The importer takes our `.gtf` genome annotation file and loads it into PostgreSQL.

It preserves:

* Chromosome / sequence ID
* Source
* Feature type
* Start position
* End position
* Score
* Strand
* Frame
* Every GTF attribute
* The original raw attribute data
* GTF metadata/header information

The database can store features such as:

```text
gene
transcript
exon
CDS
UTR
start_codon
stop_codon
...
```

We use a PostgreSQL `JSONB` column to preserve attributes that may vary between annotation types.

---

# 1. Install PostgreSQL

Download PostgreSQL from the official PostgreSQL website:

https://www.postgresql.org/download/windows/

Install PostgreSQL using the normal/default settings.

During installation, PostgreSQL will ask you to create a password for the `postgres` user.

**Remember this password.**

You will need it later.

### Important settings

Unless you intentionally changed them, use:

```text
Host:     localhost
Port:     5432
Username: postgres
```

---

# 2. Install Python

Make sure Python is installed.

Check it with:

```powershell
python --version
```

You should get something similar to:

```text
Python 3.x.x
```

If `python` doesn't work, try:

```powershell
py --version
```

---

# 3. Install the PostgreSQL Python Library

Open PowerShell in the project directory and run:

```powershell
python -m pip install "psycopg[binary]"
```

This installs `psycopg`, which allows Python to communicate with PostgreSQL.

---

# 4. Create the Database

Open **pgAdmin**.

Connect to your local PostgreSQL server.

In the left sidebar:

```text
Servers
└── PostgreSQL
    └── Databases
```

Right-click **Databases** and select:

```text
Create → Database
```

Set:

```text
Database: genome
Owner:    postgres
```

Then click **Save**.

You should now have:

```text
genome
```

as a database.

---

# 5. Get the GTF File

Make sure you have the project's NCBI GTF file.

The importer expects a file similar to:

```text
GCF_000001405.40_GRCh38.p14_genomic.gtf
```

Put the GTF somewhere accessible on your computer.

For example:

```text
project/
├── postgres_importer.py
└── data/
    └── GCF_000001405.40_GRCh38.p14_genomic.gtf
```

---

# 6. Configure the Importer

Open:

```text
postgres_importer.py
```

At the top of the file, find:

```python
DB_CONFIG = {
    "host": "localhost",
    "port": 5432,
    "dbname": "genome",
    "user": "postgres",
    "password": "YOUR_POSTGRES_PASSWORD",
}
```

Replace:

```text
YOUR_POSTGRES_PASSWORD
```

with the password you created when installing PostgreSQL.

Then find:

```python
GTF_FILE = Path(
    r"C:\path\to\your\GCF_000001405.40_GRCh38.p14_genomic.gtf"
)
```

Replace that path with the actual location of your GTF file.

For example:

```python
GTF_FILE = Path(
    r"C:\Users\Johnny\Documents\tigerhacks\data\GCF_000001405.40_GRCh38.p14_genomic.gtf"
)
```

**Important:** Keep the `r` before the path on Windows.

---

# 7. Run the Importer

Open PowerShell in the project directory.

Run:

```powershell
python postgres_importer.py
```

If everything is configured correctly, you should see something similar to:

```text
============================================================
NCBI GTF -> PostgreSQL Importer
============================================================

Database: genome
GTF:      C:\...\GCF_000001405.40_GRCh38.p14_genomic.gtf

Connected to PostgreSQL.
Creating tables...
Reading: ...
```

The importer will then load the GTF into PostgreSQL.

Depending on the size of the GTF file and the computer, this can take some time.

---

# 8. Verify the Database

After the importer finishes, it will print a summary of the imported data.

You can also open pgAdmin and navigate to:

```text
genome
└── Schemas
    └── public
        └── Tables
```

You should see:

```text
annotations
gtf_metadata
```

The important table is:

```text
annotations
```

---

# 9. Test the Database

Open the Query Tool in pgAdmin.

Run:

```sql
SELECT *
FROM annotations
WHERE feature = 'gene'
  AND gene_name = 'HBB';
```

This should return the HBB gene annotation.

For a cleaner result:

```sql
SELECT
    annotation_id,
    seqid,
    start_position,
    end_position,
    strand,
    gene_id,
    gene_name,
    attributes
FROM annotations
WHERE feature = 'gene'
  AND gene_name = 'HBB';
```

---

# 10. Find Everything Associated With HBB

To see all annotation records associated with the HBB gene:

```sql
SELECT
    feature,
    seqid,
    start_position,
    end_position,
    strand,
    gene_id,
    transcript_id,
    gene_name,
    exon_number,
    attributes
FROM annotations
WHERE gene_name = 'HBB'
ORDER BY start_position;
```

This can return things such as:

```text
gene
transcript
exon
CDS
...
```

---

# Database Structure

The database intentionally uses a generic annotation table rather than separate tables for every possible feature.

```text
annotations
│
├── annotation_id
├── seqid
├── source
├── feature
├── start_position
├── end_position
├── score
├── strand
├── frame
├── attributes_raw
├── attributes
│
├── gene_id
├── transcript_id
├── gene_name
├── exon_number
├── parent_id
├── feature_id
└── protein_id
```

### `attributes`

The `attributes` column is PostgreSQL `JSONB`.

For example, an annotation might contain information such as:

```json
{
    "gene_id": "12345",
    "gene": "HBB",
    "gene_biotype": "protein_coding",
    "description": "hemoglobin subunit beta",
    "db_xref": [
        "GeneID:3043",
        "HGNC:4827"
    ]
}
```

The exact attributes depend on the annotation.

This means we don't have to redesign the database every time the GTF contains a new attribute.

---

# Useful Queries

## Find a gene by name

```sql
SELECT *
FROM annotations
WHERE feature = 'gene'
  AND gene_name = 'HBB';
```

## Find a gene by Gene ID

```sql
SELECT *
FROM annotations
WHERE gene_id = '3043';
```

## Find all transcripts for a gene

```sql
SELECT *
FROM annotations
WHERE feature = 'transcript'
  AND gene_name = 'HBB';
```

## Find all exons for a transcript

```sql
SELECT *
FROM annotations
WHERE feature = 'exon'
  AND transcript_id = 'YOUR_TRANSCRIPT_ID'
ORDER BY start_position;
```

## See what feature types exist

```sql
SELECT
    feature,
    COUNT(*) AS count
FROM annotations
GROUP BY feature
ORDER BY count DESC;
```

## Count all annotations

```sql
SELECT COUNT(*)
FROM annotations;
```

## Find genes on a chromosome

```sql
SELECT *
FROM annotations
WHERE feature = 'gene'
  AND seqid = 'NC_000011.10';
```

---

# Troubleshooting

## `ModuleNotFoundError: No module named 'psycopg'`

Run:

```powershell
python -m pip install "psycopg[binary]"
```

Then try the importer again.

---

## `connection refused` / `WinError 10061`

PostgreSQL is probably not running.

Open:

```text
Windows Services
```

and look for the PostgreSQL service.

Start it and run the importer again.

---

## `password authentication failed`

The PostgreSQL username/password in:

```python
DB_CONFIG
```

doesn't match your PostgreSQL installation.

Make sure:

```python
"user": "postgres"
```

and that the password is correct.

---

## `database "genome" does not exist`

Open pgAdmin and create a database named:

```text
genome
```

The importer does **not** create the PostgreSQL database itself.

It creates the tables inside the database.

---

## `GTF file was not found`

Check:

```python
GTF_FILE = Path(...)
```

Make sure the path points to the actual `.gtf` file.

---

# Important: Re-running the Importer

The importer recreates the annotation tables when it starts.

That means running it again will replace the existing:

```text
annotations
gtf_metadata
```

tables and import the GTF again.

This is intentional for development because it gives everyone a clean way to rebuild the database from the source GTF.

---

# Project Architecture

The current setup is:

```text
NCBI GTF
    │
    ▼
postgres_importer.py
    │
    ▼
PostgreSQL
    │
    ├── annotations
    └── gtf_metadata
```

The application can then query PostgreSQL:

```text
Application
     │
     ▼
PostgreSQL
     │
     ▼
Genome Annotations
     │
     ├── Genes
     ├── Transcripts
     ├── Exons
     ├── CDS
     └── Other Features
```

---

# Current Limitation

The GTF contains **genome annotations and coordinates**.

It does not contain the complete `A/C/G/T` DNA sequence itself.

For example, the GTF can tell us:

```text
HBB
chromosome 11
start = ...
end   = ...
```

but we still need the corresponding GRCh38 reference genome sequence to retrieve:

```text
ATGGTGCACCTGACT...
```

The reference FASTA will therefore be a separate part of the project.

---

## Quick Setup

For someone who already has PostgreSQL installed:

```powershell
python -m pip install "psycopg[binary]"
```

Create:

```text
genome
```

in PostgreSQL.

Set the password and GTF path in:

```text
postgres_importer.py
```

Then:

```powershell
python postgres_importer.py
```

Test with:

```sql
SELECT *
FROM annotations
WHERE feature = 'gene'
  AND gene_name = 'HBB';
```

If that returns the HBB annotation, the database is working.
