"""
Core molecular data structures.

This module defines the Molecule class, which stores
basic physicochemical properties and identifiers for
small molecules retrieved from external databases.
"""


from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import Descriptors


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

    Notes
    -----
    Additional descriptors are exposed
    as computed properties using RDKit.
    """

    __slots__ = ['name', 'cid', 'molweight', 'logp',
                 'tpsa', 'molformula', 'csmiles', 'heavy_atom_count', 'hbd', 'hba', 'drug_like_lipinski', 'rdkit_mol']

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
        self.csmiles = csmiles
        self.molformula = molformula
        self.molweight = molweight
        self.logp = logp
        self.tpsa = tpsa
        self.heavy_atom_count = heavy_atom_count
        self.hba = hba
        self.hbd = hbd
        self.drug_like_lipinski = drug_like_lipinski

        if self.csmiles:
            self.rdkit_mol = Chem.MolFromSmiles(self.csmiles)
        else:
            self.rdkit_mol = None

    @property
    def morgan_fp(self):
        if self.rdkit_mol:
            morgan_gen = Chem.rdFingerprintGenerator.GetMorganGenerator(
                radius=2, fpSize=2048)
            morgan = morgan_gen.GetFingerprint(self.rdkit_mol)
            return morgan
        if self.rdkit_mol is None:
            return None

    @property
    def maccs_fp(self):
        if self.rdkit_mol:
            maccs = MACCSkeys.GenMACCSKeys(self.rdkit_mol)
            return maccs
        if self.rdkit_mol is None:
            return None

    @property
    def ring_count(self):
        if self.rdkit_mol:
            rc = Descriptors.RingCount(self.rdkit_mol)
            return rc
        if self.rdkit_mol is None:
            return None

    @property
    def aromatic_rings(self):
        if self.rdkit_mol:
            ar = Descriptors.NumAromaticRings(self.rdkit_mol)
            return ar
        if self.rdkit_mol is None:
            return None

    @property
    def rotatable_bonds(self):
        if self.rdkit_mol:
            rb = Descriptors.NumRotatableBonds(self.rdkit_mol)
            return rb
        if self.rdkit_mol is None:
            return None

    @property
    def fraction_csp3(self):
        if self.rdkit_mol:
            fcsp3 = Descriptors.FractionCSP3(self.rdkit_mol)
            return fcsp3
        if self.rdkit_mol is None:
            return None

    @property
    def formal_charge(self):
        if self.rdkit_mol:
            fc = Chem.GetFormalCharge(self.rdkit_mol)
            return fc
        if self.rdkit_mol is None:
            return None

    @property
    def num_atoms(self):
        if self.rdkit_mol:
            na = self.rdkit_mol.GetNumAtoms()
            return na
        if self.rdkit_mol is None:
            return None

    @property
    def num_bonds(self):
        if self.rdkit_mol:
            nb = self.rdkit_mol.GetNumBonds()
            return nb
        if self.rdkit_mol is None:
            return None

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
            "num_atoms": self.num_atoms,
            "heavy_atom_count": self.heavy_atom_count,
            "num_bonds": self.num_bonds,
            "hba": self.hba,
            "hbd": self.hbd,
            "ring_count": self.ring_count,
            "aromatic_rings": self.aromatic_rings,
            "rotatable_bonds": self.rotatable_bonds,
            "fraction_csp3": self.fraction_csp3,
            "formal_charge": self.formal_charge,
            "drug_like_lipinski": self.drug_like_lipinski
        }

    def __str__(self):

        return f"""\n{self.name} (CID: {self.cid}):
                \nMolecular formula: {self.molformula}
                Connectivity SMILES: {self.csmiles}
                \n
                """

    def __repr__(self):

        return f"Molecule(name='{self.name}', cid={self.cid})"
