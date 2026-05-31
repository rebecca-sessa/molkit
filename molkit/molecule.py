"""
Core molecular data structures.

This module defines the Molecule class, which stores
basic physicochemical properties and identifiers for
small molecules retrieved from external databases.
"""


class Molecule():
    """
    Representation of a small molecule.

    Parameters
    ----------
    name : str
        Common compound name.
    cid : int
        PubChem Compound ID.
    molweight : float
        Molecular weight (g/mol).
    logp : float
        Octanol/water partition coefficient.
    molformula : str
        Molecular formula.
    csmiles : str
        Canonical SMILES representation.
    tpsa : float
        Topological Polar Surface Area.
    heavy_atom_count : int
        Number of non-hydrogen atoms.
    hbd : int
        Hydrogen bond donor count.
    hba : int
        Hydrogen bond acceptor count.
    drug_like_lipinski : bool | None
        Lipinski compliance flag.
    """

    __slots__ = ['name', 'cid', 'molweight', 'logp',
                 'tpsa', 'molformula', 'csmiles', 'heavy_atom_count', 'hbd', 'hba', 'drug_like_lipinski']

    def __init__(self,
                 name: str = "",
                 cid: int = 0,
                 molweight: float = 0.0,
                 logp: float = 0.0,
                 molformula: str = "",
                 csmiles: str = "",
                 tpsa: float = 0.0,
                 heavy_atom_count: int = 0,
                 hbd: int = 0,
                 hba: int = 0,
                 drug_like_lipinski=None):

        self.name = name
        self.cid = cid
        self.molformula = molformula
        self.molweight = molweight
        self.logp = logp
        self.csmiles = csmiles
        self.tpsa = tpsa
        self.heavy_atom_count = heavy_atom_count
        self.hba = hba
        self.hbd = hbd
        self.drug_like_lipinski = drug_like_lipinski

    def to_dict(self):
        """
        Convert the molecule into a dictionary.
        Returns
        -------
        dict
            Dictionary containing molecular descriptors.
        """

        return {
            "name": self.name,
            "cid": self.cid,
            "molecular_formula": self.molformula,
            "mw": self.molweight,
            "logp": self.logp,
            "smiles": self.csmiles,
            "tpsa": self.tpsa,
            "heavy_atom_count": self.heavy_atom_count,
            "hba": self.hba,
            "hbd": self.hbd,
            "drug_like_lipinski": self.drug_like_lipinski
        }

    def __str__(self):

        return f"""\n{self.name} (CID: {self.cid}):
                \nMolecular formula: {self.molformula}
                Connectivity SMILES: {self.csmiles}
                Molecular weight: {self.molweight} g/mol
                LogP: {self.logp}
                TPSA: {self.tpsa}
                Heavy Atom Count (all atoms beside H atoms): {self.heavy_atom_count}
                H Bond Acceptor Count: {self.hba}
                H Bond Donor Count: {self.hbd}
                Drug-like: {self.drug_like_lipinski}
                \n
                """

    def __repr__(self):

        return f"Molecule({self.name} ({self.cid}))"
