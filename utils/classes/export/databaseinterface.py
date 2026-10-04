import sqlite3


class BaseDatabase:
    """ Base class for sqlite3 database interfaces """

    def __init__(self, db_path):
        self.db_path = db_path

    # database table creation

    def create_table(self, table_name, iter_columns):
        """Initialize the database and create the indicated table indexed with the indicated column """

        columns = ", ".join(list(iter_columns))
        with sqlite3.connect(self.db_path) as conn:
            conn.execute(f"""CREATE TABLE IF NOT EXISTS {table_name} ({columns});""")

    # populating the database

    def add_rows(self, table_name, iter_data):
        """Save an iterable of data tuples to the database. ncol as to correspond to the length of data."""
        list_data = list(iter_data)
        if list_data:
            placeholders = ", ".join(["?"]*len(list_data[0]))
            placeholders = f"({placeholders})"
            with sqlite3.connect(self.db_path) as conn:
                conn.executemany(f"INSERT INTO {table_name} VALUES {placeholders};", list_data)

    # getters

    def get_column_names(self, table_name):
        """ Returns the list of column names for the table """
        cols = []
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute(f"""PRAGMA table_info({table_name});""")
            for tup in cursor.fetchall():
                cols.append(tup[1])
        return cols

    def get_all_rows(self, table_name):
        """ Returns all rows in database table """
        return self._simple_fetch_query(query=f"SELECT * FROM {table_name};")

    def get_rows_pattern(self, table_name, column, value):
        """ Returns all rows in database table that fit the pattern value for the column """
        return self._simple_fetch_query(query=f"SELECT * from {table_name} WHERE {column} LIKE '%{value}%';")

    def get_rows(self, table_name, column, value):
        """ Returns all rows in database table where the column is the value """
        return self._simple_fetch_query(query=f"SELECT * from {table_name} WHERE {column} = '{value}';")

    def get_rows_custom_query(self, query):
        """ Returns all rows in database table for the provided query """
        return self._simple_fetch_query(query=query)

    def get_value(self, table_name, column, key_column, key):
        """ Returns the value for the specified column """
        return self._single_fetch_query(query=f"SELECT {column} from {table_name} WHERE {key_column} = '{key}';")

    def get_values_for_column(self, table_name, column):
        """ Returns the set of values for the specified column """
        query = f"""SELECT {column} FROM {table_name};"""
        result = self._simple_fetch_query(query=query)
        if result:
            return {tup[0] for tup in result}
        return set()

    # setters

    def update_rows(self, table_name, column, new_value, key):
        """ Replaces the row by a new row """
        update_statement = f"""UPDATE {table_name} SET {column}=? WHERE feature_id = ?;"""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute(update_statement, (new_value, key))

    def update_rows_pattern(self, table_name, column, new_value, key):
        """ Replaces the row by a new row """
        update_statement = f"""UPDATE {table_name} SET {column}=? WHERE feature_id LIKE ?;"""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute(update_statement, (new_value, f"%{key}%"))

    # deletion

    def delete_rows(self, table_name, column, value):
        """ Deletes the corresponding rows """
        query = f"DELETE FROM {table_name} WHERE {column} = '{value}';"
        self._simple_alter_query(query)

    # core query functions

    def _simple_fetch_query(self, query):
        """ Core helper function to fetch data from the database, returns a list of tuples"""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute(query)
            result = cursor.fetchall()
            return result if result else None

    def _single_fetch_query(self, query):
        """ Core helper function to fetch a single value from the database """
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute(query)
            result = cursor.fetchone()
            return result[0] if result else None

    def _simple_alter_query(self, query):
        """ Core helper function to modify data from the database """
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            cursor.execute(query)

