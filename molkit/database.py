"""
Database container for Molecule objects.

Provides utilities for storing, removing,
listing and converting molecules to pandas
DataFrames.
"""

import pandas as pd
from molkit.molecule import Molecule
from molkit.api_clients.pubchem_client import _fetch_compound


class MoleculeDatabase(list):
    """
    Collection of Molecule objects.

    Extends Python's built-in list with
    chemistry-specific helper methods.
    """

    def add_molecule(self, item):
        """
        Add a Molecule object to the database.

        Parameters
        ----------
        item : Molecule
        Molecule instance to add.
        """
        if not isinstance(item, Molecule):
            raise TypeError(f"Database only accepts Molecule class objects")
        self.append(item)

    def remove_molecule(self, item):
        if item not in self:
            raise ValueError(f"{item} not found in database")
        else:
            self.remove(item)

    def list_mols(self):
        """
        Return the names of all molecules.

        Returns
        -------
        list[str]
        """
        return [m.name for m in self]

    def to_dataframe(self):
        """
        Convert database contents to a pandas DataFrame.

        Returns
        -------
        pandas.DataFrame
        """
        return pd.DataFrame([mol.to_dict() for mol in self])

    def add_compounds(self, molname):
        """
        Retrieve compounds from PubChem
        and add them to the database.

        Parameters
        ----------
        molname : str | list[str]
        """

        _fetch_compound(molname, self)

    def __repr__(self):
        return f"MoleculeDatabase(n_molecules={len(self)})"

    def __str__(self):
        return (
            f"MoleculeDatabase "
            f"containing {len(self)} molecules")

    def get(self, name):
        for mol in self:
            if mol.name == name:
                return mol
        return None
