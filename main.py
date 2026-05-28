from molkit.molecule import Molecule
from molkit.database import MoleculeDatabase
from molkit.analysis import db_norm_similarity
from molkit.utils import export_csv, import_database
from molkit.api_clients.pubchem_client import get_compound_and_data
import requests
import csv
import numpy as np
np.set_printoptions(suppress=True, precision=4)


def main():
    print("=== Molecule Toolkit ===")
    print("1. Add molecule")
    print("2. Show database")
    print("3. Run analysis")
    print("Common abbreviations in this code: mol for molecule, db for database, fp for filepath")


if __name__ == "__main__":
    main()
