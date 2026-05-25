#HELPER FUNCTIONS
import csv
from molkit.database import MoleculeDatabase
from molkit.molecule import Molecule


def db_to_csv(database, filepath):
    if not isinstance(database, MoleculeDatabase):
            raise TypeError(f"This function only accepts Molecule class objects")
 
    with open(f"{filepath}/database.csv", "w") as db:
        writer = csv.writer(db)
        writer.writerow(["Name", "Molecular weight (g/mol)", "logP"])
        for mol in database:
            writer.writerow([mol.name, mol.molweight, mol.logP])


def csv_to_db(filepath):
    db = MoleculeDatabase()

    with open (filepath) as mol_csv:
        reader = csv.reader(mol_csv)
        next(reader)
        
        for row in reader:
            print(row)
            x = row
            mol = Molecule(
            name = x[0], 
            molweight = x[1],
            logP = x[2],
            )
            db.add_molecule(mol)  
        return db    
        
        