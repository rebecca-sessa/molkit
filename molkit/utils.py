# HELPER FUNCTIONS
import csv
from molkit.database import MoleculeDatabase
from molkit.molecule import Molecule


def export_csv(database: MoleculeDatabase, file_path_name: str):
    if not isinstance(database, MoleculeDatabase):
        raise TypeError(f"This function only accepts Molecule class objects")

    with open(f"{file_path_name}", "w") as export:
        writer = csv.writer(export)
        writer.writerow(["id", "name", "molecular_formula", "smiles",
                        "molecular_weight", "logP", "tpsa", "heavy_atom_count"])
        for molecule in database:
            writer.writerow([molecule.cid, molecule.name, molecule.molformula, molecule.csmiles,
                            molecule.molweight, molecule.logP, molecule.tpsa, molecule.HeavyAtomCount])


def import_database(filepath: str):
    db = MoleculeDatabase()

    with open(filepath) as mol_csv:
        reader = csv.reader(mol_csv)
        next(reader)

        for row in reader:
            print(row)
            x = row
            mol = Molecule(
                cid=x[0],
                name=x[1],
                molformula=x[2],
                csmiles=[3],
                molweight=[4],
                logP=[5],
                tpsa=[6],
                heavyatomcount=[7]
            )
            db.add_molecule(mol)
        return db
