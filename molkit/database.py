#BUILD DATABASE
from molkit.molecule import Molecule


class MoleculeDatabase(list):
    def add_molecule(self, item):
        if not isinstance(item, Molecule):
            raise TypeError(f"Database only accepts Molecule class objects")
        self.append(item)
        print(f"{item.name} has been added to the database!")

    def remove_molecule(self, item):
        self.remove(item)
        return f"{item.name} has been removed from the database!"
    
    def __iter__(self):
        return super().__iter__()
    
    def count(self):
        count = 0
        for m in self:
            count = count + 1
        return count
    
    def list_mols(self):
        names = []
        for m in self:
            names.append(m.name)
        return names 
