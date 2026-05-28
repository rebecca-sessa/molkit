# MolKit

Mini cheminformatics toolkit written in Python.

## Features

- Molecule objects
- Molecular database
- Molecular weight analysis
- logP filtering
- CSV import/export
- NumPy-based molecular analysis
- Feature normalization and scaling
- Outlier detection
- Pairwise similarity matrices

## Structure section

molkit/
│
├── molecule.py
├── database.py
├── analysis.py
├── utils.py
└── api_clients/

## Dependencies

- NumPy, requests

pip install -r requirements.txt


## Examples

Common abbreviations are:

- db for database
- mol(s) for molecule(s)
- fp for filepath

```python
#molecule.py
ibu = Molecule(name="Ibuprofen", molweight=200, logP=2.1)

#database.py
db1 = MoleculeDatabase()
get_compound_and_data(molname, moldb)

#analysis.py
scaled, similarity = db_normalization(mol_database)

#utils.py
export_csv(database, filepath)
