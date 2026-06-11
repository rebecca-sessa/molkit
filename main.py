"""
MolKit example workflow.

Demonstrates:

- compound retrieval from PubChem
- molecular database management
- drug-likeness analysis
- reporting
- dataset preparation
- CSV export
"""

from molkit.database import MoleculeDatabase

from molkit.analysis import calculate_lipinski
from molkit.reports import report

from molkit.datasets import prepare_ml_dataset
from molkit.io import export_csv


def main():

    db = MoleculeDatabase()

    db.add_compounds([
        "aspirin",
        "caffeine",
        "ibuprofen",
        "paracetamol"
    ])

    calculate_lipinski(db)

    print(report(db))

    print("\nDatabase preview:")
    print(db.to_dataframe().head())

    X_train, X_test, y_train, y_test = prepare_ml_dataset(db)

    print("\nTraining set shape:")
    print(X_train.shape)

    export_csv(db, "molecules.csv")


if __name__ == "__main__":
    main()
