# MolKit

MolKit is an educational cheminformatics toolkit written in Python for learning molecular data analysis, chemical informatics, visualization and machine-learning dataset preparation.

The project was developed as a hands-on learning exercise in Python, cheminformatics and computational drug discovery.

---

## Features

### Molecular representation

- Molecule objects
- Canonical SMILES support
- PubChem Compound IDs
- RDKit integration
- Morgan fingerprints
- MACCS fingerprints

### Molecular databases

- Store collections of molecules
- Add and remove compounds
- Retrieve compounds from PubChem
- Convert databases into Pandas DataFrames

### Drug-likeness analysis

- Lipinski Rule of Five
- Veber criteria
- Drug-likeness reports

### Similarity analysis

- Tanimoto similarity
- Similarity search
- Pairwise similarity matrices

### Descriptor analysis

- Descriptor correlation matrices
- Correlation ranking
- Descriptor distributions
- Pairplots
- Outlier detection
- Descriptor summaries

### Chemical space analysis

- PCA projection
- Molecular clustering
- Cluster representative selection

### Scaffold analysis

- Bemis-Murcko scaffold extraction
- Scaffold frequency analysis

### Visualization

- Similarity heatmaps
- PCA plots
- Cluster plots
- Similarity networks
- Substructure highlighting
- Murcko scaffold highlighting

### Data import/export

- CSV export
- CSV import
- Database reconstruction

### Dataset preparation

- Train/test split generation
- Feature matrix generation
- Target vector generation
- Machine-learning dataset preparation

---

## Project Structure

```text
molkit/
│
├── README.md
├── ROADMAP.md
├── pyproject.toml
├── requirements.txt
│
└── molkit/
    ├── __init__.py
    ├── molecule.py
    ├── database.py
    ├── pubchem_client.py
    ├── analysis.py
    ├── visualization.py
    ├── descriptors.py
    ├── reports.py
    ├── io.py
    └── datasets.py
```

---

## Installation

Clone the repository:

```bash
git clone https://github.com/rebecca-sessa/molkit.git
cd molkit
```

Create a conda environment:

```bash
conda create -n molkit python=3.13
conda activate molkit
```

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## Quick Start

```python
from molkit.database import MoleculeDatabase

db = MoleculeDatabase()

db.add_compounds([
    "aspirin",
    "caffeine",
    "ibuprofen"
])

print(db)
```

---

## Drug-Likeness Analysis

```python
from molkit.analysis import calculate_lipinski

calculate_lipinski(db)
```

---

## Similarity Search

```python
from molkit.analysis import similarity_search

results = similarity_search(
    db,
    target_name="aspirin",
    top_n=3
)

print(results)
```

---

## Descriptor Analysis

```python
from molkit.descriptors import descriptor_summary

summary = descriptor_summary(db)

print(summary)
```

---

## Reporting

```python
from molkit.reports import report

print(report(db))
```

---

## Dataset Preparation

```python
from molkit.analysis import calculate_lipinski
from molkit.datasets import prepare_ml_dataset

calculate_lipinski(db)

X_train, X_test, y_train, y_test = prepare_ml_dataset(db)
```

---

## Educational Purpose

MolKit is primarily intended as a learning project exploring:

- Python programming
- Cheminformatics
- Molecular descriptors
- Data analysis
- Scientific computing
- Computational drug discovery
- Machine learning workflows

The code prioritizes readability and educational value over production-level optimization.

---

## Current Version

Current stable version:

**v0.2.0**

Focus areas:

- Molecular databases
- Descriptor analysis
- Similarity analysis
- Visualization
- Dataset preparation

---

## Roadmap

Planned for v0.3.0:

- Scikit-Learn models
- Molecular classification workflows
- Model evaluation
- Cross-validation
- Feature importance analysis
- Additional machine-learning utilities