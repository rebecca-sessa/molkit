"""
MolKit main entry point.

Example workflow demonstrating:
- molecule retrieval from PubChem
- database creation
- molecular analysis
- dataset generation
- CSV export/import

This script is intended for testing and development.
"""

from molkit.database import MoleculeDatabase
from molkit.analysis import report
from molkit.utils import export_csv
from molkit.api_clients.pubchem_client import fetch_compound
from molkit.datasets import build_training_set


def main():
    """
    Demonstration workflow for MolKit.
    """

    db = MoleculeDatabase()

    compounds = [
        "aspirin",
        "caffeine",
        "ibuprofen",
        "paracetamol"
    ]

    for compound in compounds:
        fetch_compound(compound, db)

    print(report(db))

    df = db.to_dataframe()

    print("\nDatabase preview:")
    print(df.head())

    X_train, X_test, y_train, y_test = build_training_set(db)

    print("\nTraining set shape:")
    print(X_train.shape)

    export_csv(db, "molecules.csv")


if __name__ == "__main__":
    main()
