"""
Descriptor analysis utilities.

Provides functions for descriptor-based
similarity analysis, correlation analysis,
distribution visualization and outlier
detection for molecular datasets.
"""

import matplotlib.pyplot as plt
import seaborn as sns
from molkit.database import MoleculeDatabase
import pandas as pd
import numpy as np
np.set_printoptions(suppress=True, precision=4)


def descriptor_similarity_matrix(moldb: MoleculeDatabase, z_threshold=1.5):
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
            molecule.heavy_atom_count,
            molecule.ring_count,
            molecule.aromatic_rings,
            molecule.rotatable_bonds,
            molecule.fraction_csp3
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


def descriptor_correlation_matrix(database: MoleculeDatabase):
    df = database.to_dataframe()
    df_num = df.select_dtypes(include=np.number)
    df_corr = df_num.corr()

    return df_corr


def strongest_correlations(database: MoleculeDatabase, top_n=10):
    df = descriptor_correlation_matrix(database)

    mask = np.triu(np.ones(df.shape), k=1).astype(bool)
    df_flt = df.where(mask)

    df_flat = df_flt.stack()

    df_sorted = df_flat.sort_values(key=lambda x: np.abs(x), ascending=False)

    return [(desc1, desc2, valore) for (desc1, desc2), valore in df_sorted.head(top_n).items()]


def descriptor_distribution(database: MoleculeDatabase, descriptor):
    df = database.to_dataframe()

    plt.figure(figsize=(7, 7))
    sns.histplot(data=df, x=descriptor, kde=True)
    plt.title(f"{descriptor} distribution")
    plt.show()


def descriptor_pairplot(database: MoleculeDatabase):
    df = database.to_dataframe()
    df_num = df.select_dtypes(include=np.number)

    sns.pairplot(df_num)
    plt.show()


def descriptor_outliers(database, z_threshold=3):
    df = database.to_dataframe()
    fin = df.drop(columns="cid")
    df_num = fin.select_dtypes(include=np.number)

    mean_col = np.mean(df_num, axis=0)
    std_col = np.std(df_num, axis=0)

    zscore = (df_num - mean_col) / std_col

    mask = np.abs(zscore) > z_threshold

    outlier_positions = mask.stack()[mask.stack()]

    results = []
    for (idx, descriptor) in outlier_positions.index:
        results.append({
            "name": df.loc[idx, "name"],
            "descriptor": descriptor,
            "observed_value": df_num.loc[idx, descriptor],
            "zscore": zscore.loc[idx, descriptor]
        })

    return pd.DataFrame(results)


def descriptor_summary(database: MoleculeDatabase):
    df = database.to_dataframe()
    return df.describe()
