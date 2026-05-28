from molkit.database import MoleculeDatabase
from molkit.molecule import Molecule
import numpy as np
np.set_printoptions(suppress=True, precision=4)


def db_norm_similarity(moldb: MoleculeDatabase):

    data = []

    for molecule in moldb:
        data.append([
            molecule.molweight,
            molecule.logP,
            molecule.tpsa,
            molecule.HeavyAtomCount
        ])

    data = np.array(data)

    # stats
    mean_col = np.mean(data, axis=0)
    std_col = np.std(data, axis=0)

    # standardizzare z score
    zscore = (data - mean_col) / std_col

    # filtrare outlier
    mask = np.all(np.abs(zscore) <= 1.5, axis=1)
    filtered_data = data[mask]

    # normalizzare min max
    mins = filtered_data.min(axis=0)
    maxs = filtered_data.max(axis=0)

    scaled = (filtered_data - mins) / (maxs - mins)

    # matrice pairwise
    pairwise = scaled[:, None] - scaled

    # euclidean distances
    distances = np.linalg.norm(pairwise, axis=2)

    # similarity matrix
    similarity = 1 / (1 + distances)

    return scaled, similarity
