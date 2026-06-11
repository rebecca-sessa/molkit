"""
Visualization utilities for molecular datasets.

Provides plotting and rendering functions for:

- Similarity heatmaps
- Chemical space projections
- Clustering visualizations
- Similarity networks
- Substructure highlighting
- Murcko scaffold visualization

These functions generate graphical outputs from
analysis results stored in MoleculeDatabase objects.
"""

from molkit.database import MoleculeDatabase
from molkit.molecule import Molecule
from molkit.analysis import (
    similarity_matrix, pca_coordinates, cluster_molecules, tanimoto_similarity)
from molkit.descriptors import descriptor_correlation_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds import MurckoScaffold

from IPython.display import SVG, display


def plot_similarity_heatmap(database: MoleculeDatabase):
    """
    Plot pairwise molecular similarity as a heatmap.
    """
    sim_matrix = similarity_matrix(database)

    nomi = [mol.name for mol in database]

    plt.figure(figsize=(20, 20))
    sns.heatmap(sim_matrix, xticklabels=nomi, yticklabels=nomi,
                annot=True, fmt=".2f", cmap="plasma")
    plt.show()


def draw_substructure(molecule: Molecule, smarts):
    """
    Visualize a SMARTS substructure match in a molecule.
    """
    pattern = Chem.MolFromSmarts(smarts)

    if pattern is None:
        raise ValueError(f"La stringa SMARTS fornita '{smarts}' non è valida!")

    struct = molecule.rdkit_mol
    subst = struct.GetSubstructMatch(pattern)

    d = rdMolDraw2D.MolDraw2DSVG(500, 500)
    rdMolDraw2D.PrepareAndDrawMolecule(d, struct, highlightAtoms=subst)
    d.FinishDrawing()

    svg_code = d.GetDrawingText()

    display(SVG(svg_code))

    return svg_code


def plot_pca(database: MoleculeDatabase):
    """
    Plot PCA coordinates of molecular fingerprints.
    """
    coords, df_pca = pca_coordinates(database)

    sns.scatterplot(data=df_pca, x="PC1", y="PC2")
    plt.title(f"Chemical space (PCA components)")


def plot_cluster_molecules(database: MoleculeDatabase, cluster_number):
    """
    Visualize molecular clusters in PCA space.
    """
    df_pca = cluster_molecules(database, cluster_number)
    sns.scatterplot(data=df_pca, x="PC1", y="PC2",
                    hue="Cluster", palette="Set1")


def plot_similarity_network(database: MoleculeDatabase, threshold):
    """
    Build and display a molecular similarity network.
    """
    G = nx.Graph()

    for mol_a in database:
        for mol_b in database:
            if mol_a.name == mol_b.name:
                continue

            sim = tanimoto_similarity(mol_a, mol_b)

            if sim >= threshold:
                G.add_edge(mol_a.name, mol_b.name)

    plt.figure(figsize=(12, 12))
    nx.draw(G, with_labels=True, node_color="skyblue",
            edge_color="gray", node_size=800, font_size=8)
    plt.show()

    return G


def draw_murcko_scaffold(database: MoleculeDatabase, molecule_name):
    """
    Highlight the Murcko scaffold of a molecule.
    """

    for m in database:
        if m.name == molecule_name:
            structure = m.rdkit_mol

            scaffold = MurckoScaffold.GetScaffoldForMol(structure)
            scaffold_smiles = MurckoScaffold.MurckoScaffoldSmiles(
                mol=structure)

            indexes = structure.GetSubstructMatch(scaffold)

            d = rdMolDraw2D.MolDraw2DSVG(500, 500)
            rdMolDraw2D.PrepareAndDrawMolecule(
                d, structure, highlightAtoms=indexes)
            d.FinishDrawing()

            codice_svg = d.GetDrawingText()
            display(SVG(codice_svg))

            return scaffold_smiles

    print(f"Molecola '{molecule_name}' non trovata nel database.")
    return None


def descriptor_heatmap(database: MoleculeDatabase):
    df = descriptor_correlation_matrix(database)

    plt.figure(figsize=(12, 12))
    sns.heatmap(df, annot=True, fmt=".2f", cmap="plasma")
    plt.title("Descriptors correlation heatmap")
    plt.show()
