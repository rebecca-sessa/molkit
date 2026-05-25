# MolKit

Mini cheminformatics toolkit written in Python.

## Features

- Molecule objects
- Molecular database
- Molecular weight analysis
- logP filtering
- CSV import/export

## Examples

Common abbreviations are:

- db for database
- mol(s) for molecule(s)
- fp for filepath

```python
#molecule.py
from molkit.molecule import Molecule
ibu = Molecule(name="Ibuprofen", molweight=200, logP=2.1)

#database.py
from molkit.database import MoleculeDatabase
db1 = MoleculeDatabase()
db1.add_mol(ibu)

#analysis.py
from molkit.analysis import average_db, heaviest_db, lipophilic_mol
average_db = average_MW(db1)
print(average_db)

#utils.py
from molkit.utils import db_to_csv, csv_to_db
fp_db = "user/path/of/choice"
db_to_csv(db1, fp_db)
