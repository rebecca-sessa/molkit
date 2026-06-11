"""
MolKit

Educational cheminformatics toolkit for molecular
data retrieval, analysis, visualization and
machine-learning dataset preparation.
"""

from .molecule import Molecule
from .database import MoleculeDatabase

__version__ = "0.2.0"

__all__ = [
    "Molecule",
    "MoleculeDatabase"
]
