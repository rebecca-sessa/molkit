"""
Import and export utilities.

Provides CSV serialization and deserialization
for MoleculeDatabase objects.
"""

import csv
from molkit.database import MoleculeDatabase
from molkit.molecule import Molecule


def export_csv(database: MoleculeDatabase, file_path_name: str):
    """
    Export a molecular database to CSV.

    Parameters
    ----------
    database : MoleculeDatabase
    file_path_name : str
    """

    if not isinstance(database, MoleculeDatabase):

        raise TypeError(f"This function only accepts MoleculeDatabase objects")

    with open(f"{file_path_name}", "w", newline="") as export:

        writer = csv.writer(export)
        writer.writerow(["id", "name", "molecular_formula", "smiles",
                        "molecular_weight", "logP", "tpsa", "heavy_atom_count", "hba", "hbd", "drug_like_lipinski"])

        for molecule in database:
            writer.writerow([molecule.cid, molecule.name, molecule.molformula, molecule.csmiles,
                            molecule.molweight, molecule.logp, molecule.tpsa, molecule.heavy_atom_count, molecule.hba, molecule.hbd, molecule.drug_like_lipinski])


def import_database(filepath: str):
    """
    Load a molecular database from a CSV file.

    Parameters
    ----------
    filepath : str

    Returns
    -------
    MoleculeDatabase
    """
    db = MoleculeDatabase()

    with open(filepath) as mol_csv:
        reader = csv.reader(mol_csv)
        next(reader)

        for row in reader:
            x = row
            mol = Molecule(
                cid=int(x[0]),
                name=x[1],
                molformula=x[2],
                csmiles=x[3],
                molweight=float(x[4]),
                logp=float(x[5]),
                tpsa=float(x[6]),
                heavy_atom_count=int(x[7]),
                hba=int(x[8]),
                hbd=int(x[9]),
                drug_like_lipinski=(x[10].strip().lower() == "true")
            )

            db.add_molecule(mol)

        return db
