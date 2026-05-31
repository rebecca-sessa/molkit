# MolKit

MolKit is a lightweight cheminformatics toolkit written in Python for learning and experimenting with molecular data analysis, dataset preparation, and basic machine learning workflows.

This is my first Python project, developed as an educational toolkit to explore concepts in:

- Cheminformatics
- Computational drug discovery
- Molecular descriptors
- Data analysis with NumPy and Pandas
- Machine learning dataset preparation

---

## Features

### Molecular representation

- Molecule objects with common physicochemical descriptors
- Canonical SMILES support
- PubChem Compound ID integration

### Molecular databases

- Store and manage collections of molecules
- Add and remove compounds
- Convert molecular collections into Pandas DataFrames

### PubChem integration

Retrieve molecular descriptors directly from PubChem:

- Molecular weight
- Molecular formula
- XLogP
- TPSA
- Heavy atom count
- Hydrogen bond donor count
- Hydrogen bond acceptor count
- Canonical SMILES

### Molecular analysis

- Lipinski Rule of Five evaluation
- Molecular descriptor normalization
- Outlier filtering using Z-scores
- Pairwise molecular similarity matrices

### Dataset preparation

Create machine learning ready datasets:

- Training/test split
- Feature matrix generation
- Target vector generation

### Data import/export

- CSV export
- CSV import
- Database reconstruction from saved files

---

## Project structure

```text
molkit/
│
├── main.py
├── pyproject.toml
├── requirements.txt
├── README.md
│
└── molkit/
    ├── __init__.py
    ├── molecule.py
    ├── database.py
    ├── analysis.py
    ├── datasets.py
    ├── utils.py
    │
    └── api_clients/
        └── pubchem_client.py
```

---

## Installation

Clone the repository:

```bash
git clone https://github.com/rebecca-sessa/molkit.git
cd molkit
```

Create and activate a virtual environment:

```bash
conda create -n molkit python=3.13
conda activate molkit
```

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## Quick start

```python
from molkit.database import MoleculeDatabase
from molkit.api_clients.pubchem_client import fetch_compound
from molkit.analysis import report

db = MoleculeDatabase()

fetch_compound("aspirin", db)
fetch_compound("caffeine", db)

print(report(db))
```

---

## Example: DataFrame conversion

```python
df = db.to_dataframe()

print(df.head())
```

---

## Example: Machine Learning dataset

```python
from molkit.datasets import build_training_set

X_train, X_test, y_train, y_test = build_training_set(db)
```

---

## Current development status

Current version focuses on:

- Molecular descriptors
- Data analysis
- Dataset engineering

## Roadmap

Planned future developments:

- Improved molecular database reporting
- Advanced Pandas-based analysis
- Molecular descriptor visualization
- Correlation analysis between descriptors
- Machine learning workflows
- Additional molecular data sources

---

## Educational purpose

MolKit is primarily intended as a learning project to explore Python programming, cheminformatics, data science, and computational drug discovery workflows.
