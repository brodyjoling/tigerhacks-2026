# MongoDB Genomic Database Setup

This guide will walk you through setting up MongoDB, installing the required Python dependency, and importing the genomic `.gtf` data into your local MongoDB database.

## Prerequisites

Before getting started, make sure you have:

* **Python 3** installed
* The project's `mongodb_importer.py` file
* The `genomic.gtf` file
* An internet connection
* A Windows computer

---

## Step 1 — Install PyMongo

PyMongo is the Python library that allows our Python code to communicate with MongoDB.

Open a terminal in the project directory and run:

```bash
pip install pymongo
```

You should see a successful installation message when it finishes.

---

## Step 2 — Download MongoDB Community Edition

Go to the official MongoDB Community Server download page:

https://www.mongodb.com/try/download/community

Select the appropriate package for your computer.

For Windows, choose the **MSI** installer.

---

## Step 3 — Install MongoDB

Once the `.msi` file has finished downloading:

1. Open the installer.
2. Follow the installation wizard.
3. Keep the default settings unless you have a specific reason to change them.
4. Continue through the installation.
5. Wait for the installation to finish.

The installation may take several minutes.

---

## Step 4 — Open MongoDB and Create a Connection

Open **MongoDB Compass**.

Create a new connection:

1. Click **Add New Connection**.
2. Enter the connection information.
3. Give the connection a recognizable name.
4. Save/connect to it.

For a normal local MongoDB installation, the connection will typically use:

```text
mongodb://localhost:27017
```

If the project already provides a specific connection string, use that instead.

---

## Step 5 — Configure the Genomic Data Importer

Open:

```text
mongodb_importer.py
```

Near the bottom of the file, you should find the path to the genomic GTF file.

Change it so that it points to the location of your `genomic.gtf` file.

For example:

```python
gtf_path = r"C:\Users\YourName\Documents\genomic.gtf"
```

> **Important:** Using `r"..."` before the Windows path prevents backslashes from being interpreted as escape characters.

Make sure the path points to the actual `.gtf` file.

---

## Step 6 — Import the Genomic Data

Run the importer from your terminal:

```bash
python mongodb_importer.py
```

The script will read the genomic GTF data and import the relevant information into MongoDB.

Depending on the size of the GTF file and your computer, this may take a little while.

---

## Step 7 — Verify the Database

Open MongoDB Compass and check your local database.

You should see the database and collections created by the importer.

From there, you can query the genomic data and use it in the rest of the project.

---

# Done!

If everything worked correctly, you now have:

**GTF File → Python Importer → MongoDB → Your Application**

The MongoDB database can now be used to quickly look up genomic information such as genes, transcripts, exons, introns, and other annotated features without repeatedly scanning the entire GTF file.

## Troubleshooting

### `ModuleNotFoundError: No module named 'pymongo'`

Run:

```bash
pip install pymongo
```

If you have multiple Python installations, try:

```bash
python -m pip install pymongo
```

or:

```bash
py -m pip install pymongo
```

### MongoDB connection fails

Make sure MongoDB is installed and running, and verify that your connection string is correct.

For a standard local installation:

```text
mongodb://localhost:27017
```

### `FileNotFoundError`

Check the path inside `mongodb_importer.py`.

Make sure:

* The file actually exists.
* The filename is correct.
* The `.gtf` extension is included.
* Windows backslashes are handled correctly.

Example:

```python
r"C:\Users\YourName\Documents\genomic.gtf"
```

### Import is taking a long time

Large genomic GTF files contain a very large amount of annotation data. The initial import can take several minutes depending on the size of the file and your computer.

Once the data is indexed in MongoDB, lookups should be significantly faster than repeatedly scanning the entire GTF file.
