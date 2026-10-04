import sqlite3


class EdgeDatabase:
    """ Saves edges as they are produced"""

    def __init__(self, db_path):
        self.db_path = db_path
        self._create_db()

    def _create_db(self):
        """Initialize the database and create the edges table."""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute("""CREATE TABLE IF NOT EXISTS edges (node1 TEXT, node2 TEXT, edgetype TEXT)""")

    def save_edges(self, edges):
        """Save a list of (node1, node2, edgetype) tuples to the database."""
        with sqlite3.connect(self.db_path) as conn:
            conn.executemany("INSERT INTO edges (node1, node2, edgetype) VALUES (?, ?, ?)", edges)

    def load_edges(self):
        """ retrieves all edges as a list of tuples (node1, node2, edgetype)"""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute("SELECT node1, node2, edgetype FROM edges")
            return cursor.fetchall()

    def generate_edges(self, chunk_size=10000):
        """ Generates edges from chunks"""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute("SELECT node1, node2, edgetype FROM edges")
            while True:
                rows = cursor.fetchmany(chunk_size)  # Fetch in chunks
                if not rows:
                    break
                for row in rows:
                    yield row
