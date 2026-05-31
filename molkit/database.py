"""
Database container for Molecule objects.

Provides utilities for storing, removing,
listing and converting molecules to pandas
DataFrames.
"""

import pandas as pd
from molkit.molecule import Molecule


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
        try:
            self.remove(item)
            return f"{item.name} has been removed from the database!"
        except ValueError:
            return f"{item.name} is not in {self} and cannot be removed."

    def list_mols(self):
        """
        Return the names of all molecules.

        Returns
        -------
        list[str]
        """
        names = []
        for m in self:
            names.append(m.name)
        return names

    def to_dataframe(self):
        """
        Convert database contents to a pandas DataFrame.

        Returns
        -------
        pandas.DataFrame
        """
        return pd.DataFrame([mol.to_dict() for mol in self])
