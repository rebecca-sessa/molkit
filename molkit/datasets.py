"""
Dataset preparation utilities.

Functions for converting molecular databases
into machine-learning ready datasets.
"""

from molkit.database import MoleculeDatabase
import numpy as np
np.set_printoptions(suppress=True, precision=4)


def build_training_set(database: MoleculeDatabase, test_size: float = 0.2):
    """
    Split molecular descriptors into training
    and testing datasets.

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

    Returns
    -------
    tuple
        X_train, X_test, y_train, y_test
    """

    data = []

    for mol in database:
        data.append([mol.molweight,
                     mol.logp,
                     mol.hba,
                     mol.hbd,
                     mol.drug_like_lipinski])

    data_arr = np.array(data)

    n_test = int(len(data_arr) * test_size)

    rng = np.random.default_rng()

    test_indices = rng.choice(len(data_arr), size=n_test, replace=False)

    mask_test = np.zeros(len(data_arr), dtype=bool)
    mask_test[test_indices] = True

    mol_test = data_arr[mask_test]
    mol_train = data_arr[~mask_test]

    X_test = mol_test[:, :4]
    y_test = mol_test[:, 4]

    X_train = mol_train[:, :4]
    y_train = mol_train[:, 4]

    return X_train, X_test, y_train, y_test
