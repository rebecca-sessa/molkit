"""
PubChem API interface.

Functions for retrieving molecular descriptors
from the PubChem REST API and converting them
into Molecule objects.
"""
from molkit.molecule import Molecule
from molkit.database import MoleculeDatabase
import requests


def fetch_compound(molname, moldb: MoleculeDatabase):
    """
    Retrieve a compound from PubChem and add it
    to a molecular database.

    Parameters
    ----------
    molname : str
        Compound name.
    moldb : MoleculeDatabase

    Returns
    -------
    None
    """

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{molname}/property/HBondDonorCount,HBondAcceptorCount,HeavyAtomCount,TPSA,HeavyAtomCount,XLogP,MolecularWeight,MolecularFormula,CanonicalSMILES/JSON"

    headers = {
        "Accept": "application/json"
    }

    try:
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()
        result = response.json()
        compound_data = result["PropertyTable"]["Properties"][0]

        mol = Molecule(name=f"{molname}",
                       cid=int(compound_data["CID"]),
                       molweight=float(compound_data["MolecularWeight"]),
                       logp=float(compound_data.get("XLogP", 0)),
                       molformula=str(compound_data["MolecularFormula"]),
                       csmiles=str(compound_data["ConnectivitySMILES"]),
                       tpsa=float(compound_data["TPSA"]),
                       heavy_atom_count=int(compound_data["HeavyAtomCount"]),
                       hba=int(compound_data["HBondAcceptorCount"]),
                       hbd=int(compound_data["HBondDonorCount"])
                       )

        moldb.add_molecule(mol)

    except requests.exceptions.RequestException as e:
        print(f"Request error: {e}")

    except KeyError:
        print("Unexpected JSON structure")

    return None
