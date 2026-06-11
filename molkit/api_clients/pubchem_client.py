from __future__ import annotations
"""
PubChem API interface.

Functions for retrieving molecular descriptors
from the PubChem REST API and converting them
into Molecule objects.
"""
from molkit.molecule import Molecule
import requests

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from molkit.database import MoleculeDatabase


def _fetch_compound(molname, moldb: MoleculeDatabase):
    """
    Retrieve a compound from PubChem and add it
    to a molecular database.

    Parameters
    ----------
    molname : str | list[str]
    Compound name or list of compound names.

    Returns
    -------
    None
    """

    if isinstance(molname, str):
        molname = [molname]
    elif not isinstance(molname, list):
        raise TypeError("molname must be a string or a list of strings")

    for m in molname:
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{m}"
            "/property/"
            "HBondDonorCount,"
            "HBondAcceptorCount,"
            "HeavyAtomCount,"
            "TPSA,"
            "XLogP,"
            "MolecularWeight,"
            "MolecularFormula,"
            "CanonicalSMILES/JSON"
        )

        headers = {
            "Accept": "application/json"
        }

        try:
            response = requests.get(url, headers=headers, timeout=10)
            response.raise_for_status()
            result = response.json()
            compound_data = result["PropertyTable"]["Properties"][0]

            mol = Molecule(name=m,
                           cid=int(compound_data["CID"]),
                           molweight=float(compound_data["MolecularWeight"]),
                           logp=float(compound_data.get("XLogP", 0)),
                           molformula=str(compound_data["MolecularFormula"]),
                           csmiles=str(compound_data["ConnectivitySMILES"]),
                           tpsa=float(compound_data["TPSA"]),
                           heavy_atom_count=int(
                               compound_data["HeavyAtomCount"]),
                           hba=int(compound_data["HBondAcceptorCount"]),
                           hbd=int(compound_data["HBondDonorCount"])
                           )
            if m in moldb.list_mols():
                continue
            else:
                moldb.add_molecule(mol)

        except requests.exceptions.RequestException as e:
            print(f"Failed to load {m}: {e}")

        except KeyError:
            print(f"Unexpected JSON structure for {m}")
