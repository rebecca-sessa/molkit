"""

Reporting utilities.

Provides textual and tabular summaries

for molecular databases, including

drug-likeness statistics and database

overview reports.

"""
from molkit.database import MoleculeDatabase
from molkit.analysis import calculate_lipinski
import pandas as pd


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


def druglikeness_report(database: MoleculeDatabase):
    """
    Generate a combined Lipinski and Veber report.

    Parameters
    ----------
    database : MoleculeDatabase

    Returns
    -------
    pandas.DataFrame
    Drug-likeness summary for all molecules.
    """
    x = []

    for m in database:
        name = m.name
        if m.drug_like_lipinski == True:
            lip = True
        else:
            lip = False

        if m.tpsa <= 140 and m.rotatable_bonds <= 10:
            vb = True
        else:
            vb = False

        mol = {
            "name molecule": name,
            "Lipinski pass": lip,
            "Veber pass": vb,
        }

        x.append(mol)

    df = pd.DataFrame(x)
    df["Veber and Lipinski pass"] = df["Lipinski pass"] & df["Veber pass"]

    return df
