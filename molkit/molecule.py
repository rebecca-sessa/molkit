# BUILD MOLECULE

class Molecule():
    __slots__ = ['name', 'cid', 'molweight', 'logP',
                 'tpsa', 'molformula', 'csmiles', 'HeavyAtomCount']

    def __init__(self, name="", cid=0, molweight=0., logP=0., molformula="", csmiles="", tpsa=0., heavyatomcount=0.):
        self.name = name
        self.cid = cid
        self.molformula = molformula
        self.molweight = molweight
        self.logP = logP
        self.csmiles = csmiles
        self.tpsa = tpsa
        self.HeavyAtomCount = heavyatomcount

    def is_drug_like(self):
        if self.molweight <= 500:
            return True
        else:
            return False

    def __str__(self):
        return f"""\n{self.name} (CID: {self.cid}):
\nMolecular formula: {self.molformula}
Connectivity SMILES: {self.csmiles}
Molecular weight: {self.molweight} g/mol
LogP: {self.logP}
TPSA: {self.tpsa}
Heavy Atom Count (all atoms beside H atoms): {self.HeavyAtomCount}
Drug-like: {self.is_drug_like()}
\n
"""

    def __repr__(self):
        return f"Molecule({self.name} ({self.cid}))"
