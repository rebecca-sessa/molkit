# API REQUEST
from molkit.molecule import Molecule
from molkit.database import MoleculeDatabase
import requests


def get_compound_and_data(molname, moldb: MoleculeDatabase):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{molname}/property/HeavyAtomCount,TPSA,HeavyAtomCount,XLogP,MolecularWeight,MolecularFormula,CanonicalSMILES/JSON"

    headers = {
        "Accept": "application/json"
    }

    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()
        result = response.json()
        compound_data = result["PropertyTable"]["Properties"][0]

        mol = Molecule(name=f"{molname}",
                       cid=int(compound_data["CID"]),
                       molweight=float(compound_data["MolecularWeight"]),
                       logP=float(compound_data["XLogP"]),
                       molformula=str(compound_data["MolecularFormula"]),
                       csmiles=str(compound_data["ConnectivitySMILES"]),
                       tpsa=float(compound_data["TPSA"]),
                       heavyatomcount=float(compound_data["HeavyAtomCount"])
                       )

        moldb.add_molecule(mol)

        return moldb

    except requests.exceptions.RequestException as e:
        print(f"Request error: {e}")

    except KeyError:
        print("Unexpected JSON structure")

    return None
