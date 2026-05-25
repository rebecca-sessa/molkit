#ANALYSIS
from molkit.molecule import Molecule
from molkit.database import MoleculeDatabase
from datetime import datetime


def average_MW(database):
    if not isinstance(database, MoleculeDatabase):
            raise TypeError(f"This function only accepts MoleculeDatabase class objects")
    
    if len(database) == 0:
                raise ValueError("Database is empty!")

    else:
        nom = 0
        dnom = 0
        
        for m in database:
            nom = nom + m.molweight
            dnom = dnom + 1
            
        else:
            mean = nom / dnom
            dt = datetime.now() #current
        
            return f"""The average molecular weight of the database is: {mean} g/mol.
Analysis performed on: {dt}"""


def heaviest_mol(database):
    heaviest = 0
    name = ""
    
    if not isinstance(database, MoleculeDatabase):
            raise TypeError(f"This function only accepts MoleculeDatabase class objects")
    
    else:
        for m in database:
            if m.molweight > heaviest:
                heaviest = m.molweight
                name = m.name
        dt = datetime.now() #current
        
        return f"""The heaviest molecule is {name} with a molecular weight of {heaviest} g/mol.
Analysis performed on: {dt}"""


def filterby_logP(database, treshold):
    filtered = MoleculeDatabase()
    
    if not isinstance(database, MoleculeDatabase):
            raise TypeError(f"This function only accepts MoleculeDatabase class objects")
    
    else:
        for m in database:
            if m.logP > treshold:
                filtered.add_molecule(m)
        dt = datetime.now() #current
        
        return f"""Filtered molecules in {database} following {treshold} treshold:
        {filtered}
Analysis performed on: {dt}"""