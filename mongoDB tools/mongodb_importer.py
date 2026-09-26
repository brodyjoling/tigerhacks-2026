from pymongo import MongoClient

client = MongoClient("mongodb://localhost:27017/")
db = client["genome"]

collections = {
    "gene": db["genes"],
    "transcript": db["transcripts"],
    "exon": db["exons"],
    "CDS": db["cds"]
}


def parse_attributes(attribute_string):

    attributes = {}

    for item in attribute_string.strip().split(";"):

        item = item.strip()

        if not item:
            continue

        key, value = item.split(" ", 1)

        attributes[key] = value.strip('"')

    return attributes


def import_gtf(filename):

    batches = {
        "gene": [],
        "transcript": [],
        "exon": [],
        "CDS": []
    }

    BATCH_SIZE = 1000

    with open(filename, "r") as file:

        for line in file:

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            feature = fields[2]

            if feature not in batches:
                continue

            attributes = parse_attributes(fields[8])

            document = {
                "chromosome": fields[0],
                "source": fields[1],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "score": fields[5],
                "strand": fields[6],
                "frame": fields[7],
                **attributes
            }

            batches[feature].append(document)

            if len(batches[feature]) >= BATCH_SIZE:

                collections[feature].insert_many(
                    batches[feature]
                )

                batches[feature].clear()

    # Insert remaining documents
    for feature, batch in batches.items():

        if batch:
            collections[feature].insert_many(batch)

# use only once
import_gtf(r"data/ncbi_dataset/data/GCF_000001405.40/genomic.gtf")