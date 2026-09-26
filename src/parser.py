

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(PROJECT_ROOT))

from src.ProtoDNASequencer import ProtoDNASequencer


GTF_FILE = PROJECT_ROOT / "genomic.gtf"
FASTA_FILE = PROJECT_ROOT / "GCF_000001405.40_GRCh38.p14_genomic.fna"

sequencer = ProtoDNASequencer(FASTA_FILE)
sequencer.load_sequence()


def parse_attributes(attribute_string):
    attributes = {}

    for item in attribute_string.split(";"):
        item = item.strip()

        if not item:
            continue

        parts = item.split(" ", 1)

        if len(parts) == 1:
            attributes[parts[0]] = True
            continue

        key, value = parts
        attributes[key] = value.strip('"')

    return attributes
    
def find_gene_by_name(gene_name):

    with open(GTF_FILE, "r") as file:

        for line in file:

            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")

            if fields[2] != "gene":
                continue

            attributes = parse_attributes(fields[8])

            if attributes.get("gene_name", attributes.get("gene")) == gene_name:
                return {
                    "chromosome": fields[0],
                    "gene": attributes.get("gene_name", attributes.get("gene")),
                    "gene_id": attributes.get("gene_id"),
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6]
                }

    return None

def find_gene(chromosome, position):

    with open(GTF_FILE, "r") as file:

        for line in file:

            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")

            if fields[0] != chromosome:
                continue

            if fields[2] != "gene":
                continue

            start = int(fields[3])
            end = int(fields[4])

            if start <= position <= end:

                attributes = parse_attributes(fields[8])

                return {
                    "gene": attributes.get("gene_name", attributes.get("gene")),
                    "chromosome": fields[0],
                    "gene_id": attributes.get("gene_id"),
                    "start": start,
                    "end": end,
                    "strand": fields[6]
                }

    return None
def find_transcripts(gene_id):

    transcripts = []

    with open(GTF_FILE, "r") as file:

        for line in file:

            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")

            if fields[2] != "transcript":
                continue

            attributes = parse_attributes(fields[8])

            if attributes.get("gene_id") != gene_id:
                continue

            transcripts.append({
                "transcript_id": attributes.get("transcript_id"),
                "transcript_name": attributes.get("transcript_name"),
                "gene_id": gene_id,
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6]
            })

    return transcripts
def find_exons(transcript_id):

    exons = []

    with open(GTF_FILE, "r") as file:

        for line in file:

            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")

            if fields[2] != "exon":
                continue

            attributes = parse_attributes(fields[8])

            if attributes.get("transcript_id") != transcript_id:
                continue

            exons.append({
                "exon_id": attributes.get("exon_id"),
                "transcript_id": transcript_id,
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6]
            })

    return exons

def get_sequence(chromosome, start, end):
    first_chromosome = sequencer.getRecord().id

    if chromosome != first_chromosome:
        return None

    # GTF coordinates are 1-based and inclusive.
    # Python slices are 0-based and end-exclusive.
    return sequencer.getSeq()[start - 1:end]

first_chromosome_gene = find_gene("NC_000001.11", 14000)

print(first_chromosome_gene)

if first_chromosome_gene:
    print(get_sequence(
        first_chromosome_gene["chromosome"],
        first_chromosome_gene["start"],
        first_chromosome_gene["end"],
    ))
# gene = find_gene("chr1", 14001)
# print(gene)
# print("\n")

# transcripts = find_transcripts(gene["gene_id"])
# print(transcripts)
# print("\n")

# exons = find_exons(transcripts[0]["transcript_id"])
# print(exons)
# print("\n")

# print(get_sequence("chr1", exons[0]["start"], exons[0]["end"]))

# [{'transcript_id': 'ENST00000832824.1', 'transcript_name': 'DDX11L16-260', 'gene_id': 'ENSG00000290825.2', 'start': 11121, 'end': 14413, 'strand': '+'},
#   {'transcript_id': 'ENST00000832825.1', 'transcript_name': 'DDX11L16-261', 'gene_id': 'ENSG00000290825.2', 'start': 11125, 'end': 14405, 'strand': '+'},
#     {'transcript_id': 'ENST00000832826.1', 'transcript_name': 'DDX11L16-262', 'gene_id': 'ENSG00000290825.2', 'start': 11410, 'end': 14413, 'strand': '+'},
#       {'transcript_id': 'ENST00000832827.1', 'transcript_name': 'DDX11L16-263', 'gene_id': 'ENSG00000290825.2', 'start': 11411, 'end': 14413, 'strand': '+'},
#         {'transcript_id': 'ENST00000832828.1', 'transcript_name': 'DDX11L16-264', 'gene_id': 'ENSG00000290825.2', 'start': 11426, 'end': 14409, 'strand': '+'},
#           {'transcript_id': 'ENST00000832829.1', 'transcript_name': 'DDX11L16-265', 'gene_id': 'ENSG00000290825.2', 'start': 11770, 'end': 14416, 'strand': '+'},
#             {'transcript_id': 'ENST00000832830.1', 'transcript_name': 'DDX11L16-266', 'gene_id': 'ENSG00000290825.2', 'start': 11819, 'end': 14413, 'strand': '+'},
#               {'transcript_id': 'ENST00000832837.1', 'transcript_name': 'DDX11L16-273', 'gene_id': 'ENSG00000290825.2', 'start': 11823, 'end': 14406, 'strand': '+'},
#                 {'transcript_id': 'ENST00000832836.1', 'transcript_name': 'DDX11L16-272', 'gene_id': 'ENSG00000290825.2', 'start': 11824, 'end': 14409, 'strand': '+'},
#                   {'transcript_id': 'ENST00000832832.1', 'transcript_name': 'DDX11L16-268', 'gene_id': 'ENSG00000290825.2', 'start': 11824, 'end': 14413, 'strand': '+'},
#                     {'transcript_id': 'ENST00000832833.1', 'transcript_name': 'DDX11L16-269', 'gene_id': 'ENSG00000290825.2', 'start': 11824, 'end': 14413, 'strand': '+'}, 
#                     {'transcript_id': 'ENST00000832831.1', 'transcript_name': 'DDX11L16-267', 'gene_id': 'ENSG00000290825.2', 'start': 11824, 'end': 14416, 'strand': '+'},
#                       {'transcript_id': 'ENST00000832834.1', 'transcript_name': 'DDX11L16-270', 'gene_id': 'ENSG00000290825.2', 'start': 11825, 'end': 14413, 'strand': '+'}, 
#                       {'transcript_id': 'ENST00000832835.1', 'transcript_name': 'DDX11L16-271', 'gene_id': 'ENSG00000290825.2', 'start': 11828, 'end': 14416, 'strand': '+'}, 
#                       {'transcript_id': 'ENST00000832839.1', 'transcript_name': 'DDX11L16-275', 'gene_id': 'ENSG00000290825.2', 'start': 11845, 'end': 14417, 'strand': '+'}, 
#                       {'transcript_id': 'ENST00000832841.1', 'transcript_name': 'DDX11L16-277', 'gene_id': 'ENSG00000290825.2', 'start': 11847, 'end': 14415, 'strand': '+'}, 
#                       {'transcript_id': 'ENST00000832840.1', 'transcript_name': 'DDX11L16-276', 'gene_id': 'ENSG00000290825.2', 'start': 11847, 'end': 14416, 'strand': '+'}, 
#                       {'transcript_id': 'ENST00000832838.1', 'transcript_name': 'DDX11L16-274', 'gene_id': 'ENSG00000290825.2', 'start': 11847, 'end': 14421, 'strand': '+'}, 
#                       {'transcript_id': 'ENST00000832842.1', 'transcript_name': 'DDX11L16-278', 'gene_id': 'ENSG00000290825.2', 'start': 11850, 'end': 14410, 'strand': '+'}, 
#                       {'transcript_id': 'ENST00000456328.3', 'transcript_name': 'DDX11L16-258', 'gene_id': 'ENSG00000290825.2', 'start': 11850, 'end': 14416, 'strand': '+'}, 
#                       {'transcript_id': 'ENST00000832845.1', 'transcript_name': 'DDX11L16-281', 'gene_id': 'ENSG00000290825.2', 'start': 11854, 'end': 14410, 'strand': '+'}, 
#                       {'transcript_id': 'ENST00000832843.1', 'transcript_name': 'DDX11L16-279', 'gene_id': 'ENSG00000290825.2', 'start': 11854, 'end': 14413, 'strand': '+'}, 
#                       {'transcript_id': 'ENST00000832844.1', 'transcript_name': 'DDX11L16-280', 'gene_id': 'ENSG00000290825.2', 'start': 11854, 'end': 14413, 'strand': '+'}, 
#                       {'transcript_id': 'ENST00000832847.1', 'transcript_name': 'DDX11L16-283', 'gene_id': 'ENSG00000290825.2', 'start': 11883, 'end': 14413, 'strand': '+'}, 
#                       {'transcript_id': 'ENST00000832846.1', 'transcript_name': 'DDX11L16-282', 'gene_id': 'ENSG00000290825.2', 'start': 11883, 'end': 14414, 'strand': '+'}, 
#                       {'transcript_id': 'ENST00000832848.1', 'transcript_name': 'DDX11L16-284', 'gene_id': 'ENSG00000290825.2', 'start': 12259, 'end': 14407, 'strand': '+'}, 
#                       {'transcript_id': 'ENST00000832849.1', 'transcript_name': 'DDX11L16-285', 'gene_id': 'ENSG00000290825.2', 'start': 12524, 'end': 14410, 'strand': '+'}, 
#                       {'transcript_id': 'ENST00000832823.1', 'transcript_name': 'DDX11L16-259', 'gene_id': 'ENSG00000290825.2', 'start': 14404, 'end': 24894, 'strand': '+'}]

# [
# {'exon_id': 'ENSE00004248723.1', 'transcript_id': 'ENST00000832824.1', 'start': 11121, 'end': 11211, 'strand': '+'},
# {'exon_id': 'ENSE00004248735.1', 'transcript_id': 'ENST00000832824.1', 'start': 12010, 'end': 12227, 'strand': '+'}, 
# {'exon_id': 'ENSE00003582793.1', 'transcript_id': 'ENST00000832824.1', 'start': 12613, 'end': 12721, 'strand': '+'}, 
# {'exon_id': 'ENSE00004248730.1', 'transcript_id': 'ENST00000832824.1', 'start': 13453, 'end': 14413, 'strand': '+'}
# ]
