"""
Molecular analysis utilities.

Contains functions for similarity analysis,
drug-likeness evaluation and database reporting.
"""

from molkit.database import MoleculeDatabase
import numpy as np
np.set_printoptions(suppress=True, precision=4)


def db_norm_similarity(moldb: MoleculeDatabase, z_threshold=1.5):
    """
    Compute a normalized similarity matrix.

    Molecular descriptors are standardized,
    filtered for outliers and min-max scaled
    before pairwise Euclidean distances are
    converted into similarity scores.

    Parameters
    ----------
    moldb : MoleculeDatabase
    z_threshold : float

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Scaled descriptor matrix and similarity matrix.
    """
    data = []

    for molecule in moldb:
        data.append([
            molecule.molweight,
            molecule.logp,
            molecule.tpsa,
            molecule.heavy_atom_count
        ])

    data = np.array(data)

    # stats
    mean_col = np.mean(data, axis=0)
    std_col = np.std(data, axis=0)

    # standardizzare z score
    zscore = (data - mean_col) / std_col

    # filtrare outlier
    mask = np.all(np.abs(zscore) <= z_threshold, axis=1)
    filtered_data = data[mask]

    # normalizzare min max
    mins = filtered_data.min(axis=0)
    maxs = filtered_data.max(axis=0)

    ranges = maxs - mins
    ranges[ranges == 0] = 1

    scaled = (filtered_data - mins) / ranges

    # matrice pairwise
    pairwise = scaled[:, None] - scaled

    # euclidean distances
    distances = np.linalg.norm(pairwise, axis=2)

    # similarity matrix
    similarity = 1 / (1 + distances)

    return scaled, similarity


def calculate_lipinski(database: MoleculeDatabase):
    """
    Annotate molecules according to Lipinski's Rule of Five.

    Parameters
    ----------
    database : MoleculeDatabase

    Returns
    -------
    MoleculeDatabase
    """
    for mol in database:
        if mol.molweight <= 500 and mol.logp <= 5 and mol.hbd <= 5 and mol.hba <= 10:
            mol.drug_like_lipinski = True
        else:
            mol.drug_like_lipinski = False

    return database


def report(database: MoleculeDatabase):
    """
    Generate a summary report for a molecular database.

    Returns
    -------
    str
        Human-readable report.
    """

    calculate_lipinski(database)
    df = database.to_dataframe()

    # mol count
    counter = len(df)

    # mean mw and logp
    mean_mw = df["mw"].mean()
    mean_logp = df["logp"].mean()

    # druglike
    dl = len(df[df["drug_like_lipinski"]])
    ndl = len(df[~df["drug_like_lipinski"]])

    return f"""\nReport:
    \nCount: {counter} molecules
    Average MW: {mean_mw:.3f} g/mol
    Average LogP: {mean_logp:.3f}
    Drug like: {dl} molecules
    \nNon drug like: {ndl} molecules
    """
