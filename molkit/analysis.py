"""
Molecular analysis utilities.

Contains functions for similarity analysis,
drug-likeness evaluation and database reporting.
"""
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import pandas as pd
from molkit.molecule import Molecule
from molkit.database import MoleculeDatabase
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import DataStructs
from rdkit.DataStructs import TanimotoSimilarity

import numpy as np
np.set_printoptions(suppress=True, precision=4)


"""
Molecular analysis utilities.

Provides functions for:

- Drug-likeness evaluation (Lipinski, Veber)
- Molecular similarity analysis
- Substructure searching
- Chemical space projection (PCA)
- Molecular clustering
- Scaffold analysis

These functions operate on MoleculeDatabase objects
and return numerical results, tables or molecular
annotations suitable for further processing.
"""


# Drug - likeness filtering

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


def calculate_veber(database: MoleculeDatabase):
    """
    Evaluate Veber's oral bioavailability criteria.

    Parameters
    ----------
    database : MoleculeDatabase

    Returns
    -------
    pandas.DataFrame
    DataFrame containing TPSA, rotatable bonds
    and Veber compliance for each molecule.
    """
    veber_flt = []

    for m in database:
        if m.tpsa < 140 and m.rotatable_bonds < 10:
            veber = True
        else:
            veber = False

        mol = {
            "name": m.name,
            "cid": m.cid,
            "molecular_formula": m.molformula,
            "mw": m.molweight,
            "tpsa": m.tpsa,
            "rotatable_bonds": m.rotatable_bonds,
            "rdkit object": m.rdkit_mol,
            "veber filter": veber
        }

        veber_flt.append(mol)

    df = pd.DataFrame(veber_flt)

    return df


# Similarity

def tanimoto_similarity(mol_a: Molecule, mol_b: Molecule):
    """
    Compute Tanimoto similarity between two molecules.

    Parameters
    ----------
    mol_a : Molecule
    mol_b : Molecule

    Returns
    -------
    float
    Tanimoto coefficient between Morgan fingerprints.
    """
    r = TanimotoSimilarity(mol_a.morgan_fp, mol_b.morgan_fp)
    return r


def similarity_search(database: MoleculeDatabase, target_name: str, top_n):
    """
    Find the most similar molecules to a target compound.

    Parameters
    ----------
    database : MoleculeDatabase
    target_name : str
    top_n : int

    Returns
    -------
    list[dict]
    Molecules ranked by Tanimoto similarity.
    """
    target = None

    for i in database:
        if i.name == target_name:
            target = i
            break

    if target is None:
        raise ValueError(

            f"Molecule '{target_name}' not found"

        )

    sim = []
    for i in database:
        if i is not target:
            name_dict = i.name
            sim_dict = tanimoto_similarity(i, target)
            d = {
                "name": name_dict,
                "similarity": sim_dict
            }
            sim.append(d)

    sim_sorted = sorted(sim, key=lambda x: x['similarity'], reverse=True)

    out = sim_sorted[:top_n]

    return out


def similarity_matrix(database: MoleculeDatabase):
    """
    Compute pairwise Tanimoto similarity matrix.

    Parameters
    ----------
    database : MoleculeDatabase

    Returns
    -------
    numpy.ndarray
    Square similarity matrix.
    """
    n_mol = len(database)

    m = np.zeros((n_mol, n_mol))

    for i in range(n_mol):
        for j in range(n_mol):
            mol1 = database[i]
            mol2 = database[j]
            ris = tanimoto_similarity(mol1, mol2)
            m[i, j] = ris

    return m


#  Substructure search

def substructure_search(database: MoleculeDatabase, smarts):
    """
    Search molecules containing a SMARTS pattern.

    Parameters
    ----------
    database : MoleculeDatabase
    smarts : str

    Returns
    -------
    list[Molecule]
    Molecules matching the query pattern.
    """
    res = []
    pattern = Chem.MolFromSmarts(smarts)

    for i in database:
        struct = i.rdkit_mol
        subst = struct.HasSubstructMatch(pattern)
        if subst == True:
            res.append(i)

    return res


# Dimensionality reduction

def pca_coordinates(database: MoleculeDatabase):
    """
    Project molecular fingerprints into a lower-dimensional space.

    Parameters
    ----------
    database : MoleculeDatabase

    Returns
    -------
    tuple
    PCA coordinates and corresponding DataFrame.
    """
    matrix = []

    for m in database:
        fp = m.morgan_fp
        arr = np.zeros((2048,), dtype=int)
        DataStructs.ConvertToNumpyArray(fp, arr)
        matrix.append(arr)

    X = np.array(matrix)

    pca = PCA(n_components=2)
    coords = pca.fit_transform(X)

    x_ax = coords[:, 0]
    y_ax = coords[:, 1]

    names = [m.name for m in database]

    graph = {
        "mol_name": names,
        "PC1": x_ax,
        "PC2": y_ax
    }

    df_pca = pd.DataFrame(graph)

    return coords, df_pca


# Clustering

def cluster_molecules(database: MoleculeDatabase, cluster_number):
    """
    Cluster molecules using K-Means on PCA coordinates.

    Parameters
    ----------
    database : MoleculeDatabase
    cluster_number : int

    Returns
    -------
    pandas.DataFrame
    PCA coordinates with assigned cluster labels.
    """
    coords, df_pca = pca_coordinates(database)

    kmeans = KMeans(n_clusters=cluster_number, random_state=42)
    groups = kmeans.fit_predict(coords)
    df_pca["Cluster"] = groups

    return df_pca


def cluster_representatives(database: MoleculeDatabase, cluster_number):
    """
    Identify representative molecules for each cluster.

    Parameters
    ----------
    database : MoleculeDatabase
    cluster_number : int

    Returns
    -------
    list[Molecule]
    Representative molecule of each cluster.
    """
    df_pca = cluster_molecules(database, cluster_number)
    rep = []

    for i in range(cluster_number):
        cluster = df_pca[df_pca["Cluster"] == i]
        c_names = cluster["mol_name"].tolist()
        c_mols = [m for m in database if m.name in c_names]
        best_score = -1
        represent = None

        for mol_a in c_mols:
            t = 0
            for mol_b in c_mols:
                t += tanimoto_similarity(mol_a, mol_b)

            if t > best_score:
                best_score = t
                represent = mol_a

        rep.append(represent)

    return rep


# Scaffolds

def get_murcko_scaffold(database: MoleculeDatabase, molecule_name):
    """
    Extract the Bemis-Murcko scaffold of a molecule.

    Parameters
    ----------
    database : MoleculeDatabase
    molecule_name : str

    Returns
    -------
    str | None
    Scaffold SMILES if found.
    """
    for m in database:
        if m.name == molecule_name:
            structure = m.rdkit_mol

            scaffold_smiles = MurckoScaffold.MurckoScaffoldSmiles(
                mol=structure)

            return scaffold_smiles

    print(f"Molecola '{molecule_name}' non trovata nel database.")
    return None


def scaffold_frequency(database: MoleculeDatabase):
    """
    Count scaffold occurrences in a molecular database.

    Parameters
    ----------
    database : MoleculeDatabase

    Returns
    -------
    pandas.Series
    Scaffold frequencies sorted by occurrence.
    """
    smiles = []

    for m in database:
        structure = m.rdkit_mol
        scaffold_smiles = MurckoScaffold.MurckoScaffoldSmiles(mol=structure)

        smiles.append(scaffold_smiles)

    smi = pd.Series(smiles)

    return smi.value_counts()
