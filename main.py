from molkit.molecule import Molecule
from molkit.database import MoleculeDatabase
from molkit.analysis import average_MW, heaviest_mol, filterby_logP
from molkit.utils import db_to_csv, csv_to_db
import csv

def main():
    print("=== Molecule Toolkit ===")
    print("1. Add molecule")
    print("2. Show database")
    print("3. Run analysis")
    print("Common abbreviations in this code: mol for molecule, db for database, fp for filepath")

if __name__ == "__main__":
    main()

