import psycopg
from psycopg.rows import dict_row

DB_CONFIG = {
    "host": "localhost",
    "port": 5432,
    "dbname": "genome",
    "user": "postgres",
    "password": "password",
}

with psycopg.connect(**DB_CONFIG, row_factory=dict_row) as connection:
    with connection.cursor() as cursor:
        cursor.execute("SELECT * FROM annotations LIMIT 1;")
        row = cursor.fetchone()
        print(row)