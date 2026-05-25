#BUILD MOLECULE

class Molecule():
    __slots__ = ['name', 'molweight', 'logP']
    
    def __init__(self, name="default", molweight=0, logP=0):
        self.name = name
        self.molweight = molweight
        self.logP = logP
    
    def is_drug_like(self):
        if self.molweight <= 500:
            return True
        else:
            return False
    
    def __str__(self):
        return f"""
{self.name}:
Molecular weight: {self.molweight} g/mol
LogP: {self.logP}
Drug-like: {self.is_drug_like()}
"""
    
    def __repr__(self):
        return f"Molecule({self.name}, {self.molweight}, {self.logP})"
    


