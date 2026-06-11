"""
Dataset preparation utilities.

Functions for converting molecular databases
into machine-learning ready datasets.
"""

from molkit.database import MoleculeDatabase
import numpy as np
np.set_printoptions(suppress=True, precision=4)


def prepare_ml_dataset(database: MoleculeDatabase, test_size: float = 0.2, random_state: int = 42):
    """
    Prepare a machine-learning dataset from a molecular database.

    Features
    --------
    - molecular weight
    - logP
    - HBA
    - HBD

    Target
    ------
    - drug_like_lipinski

    Parameters
    ----------
    database : MoleculeDatabase
    test_size : float
    random_state : int

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        X_train, X_test, y_train, y_test
    """

    if any(m.drug_like_lipinski is None for m in database):
        raise ValueError(
            "Run calculate_lipinski() before preparing the dataset."
        )

    data = []

    for mol in database:
        data.append([
            mol.molweight,
            mol.logp,
            mol.hba,
            mol.hbd,
            int(mol.drug_like_lipinski)
        ])

    data = np.array(data, dtype=float)

    n_test = int(len(data) * test_size)

    rng = np.random.default_rng(random_state)

    test_indices = rng.choice(
        len(data),
        size=n_test,
        replace=False
    )

    mask = np.zeros(len(data), dtype=bool)
    mask[test_indices] = True

    train = data[~mask]
    test = data[mask]

    X_train = train[:, :4]
    y_train = train[:, 4].astype(int)

    X_test = test[:, :4]
    y_test = test[:, 4].astype(int)

    return X_train, X_test, y_train, y_test
