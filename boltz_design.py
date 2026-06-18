import os
from boltz_api import Boltz

boltz_key = os.getenv('BOLTZ_API_KEY')
print('got key')
client = Boltz(api_key=boltz_key)
print('used key')

target = {
    "entities": [
      # protein chains only; at least one
      { "type": "protein", "value": "MKTIIALSYIFCLVFA", "chain_ids": ["A"] }
    ],
    "pocket_residues": { "A": [2, 3, 4, 7, 8, 9] }, # optional: keyed by chain ID; pocket residues (0-indexed); omit to auto-detect
    "reference_ligands": ["CC(=O)Oc1ccccc1C(=O)O"], # optional: known-binder SMILES that help locate the pocket
    "constraints": [ # optional: guide the geometry
      {
        "type": "pocket", # keep the binder near a set of receptor residues
        "binder_chain_id": "L", # chain ID the pipeline assigns to the designed molecule
        "contact_residues": { "A": [2, 3, 4, 7, 8, 9] },
        "max_distance_angstrom": 6.0
      }
    ]
  }

molecule_filters =  {
    "boltz_smarts_catalog_filter_level": "recommended", # recommended | extra | aggressive | disabled
    "custom_filters": [
      # any combination; a molecule must pass all of them (AND logic). one of each type shown:
      {
        "type": "lipinski_filter", # Rule of Five
        "max_mw": 500,
        "max_logp": 5,
        "max_hbd": 5,
        "max_hba": 10,
        "allow_single_violation": False # optional: allow one rule to fail
      },
      {
        "type": "rdkit_descriptor_filter", # min/max on RDKit descriptors; include only the ones you want to bound
        "mol_wt": { "min": 150, "max": 500 },
        "mol_logp": { "max": 5 },
        "tpsa": { "max": 140 },
        "num_h_donors": { "max": 5 },
        "num_h_acceptors": { "max": 10 },
        "num_rotatable_bonds": { "max": 10 },
        "num_heteroatoms": { "max": 12 },
        "num_aromatic_rings": { "min": 1, "max": 4 },
        "num_rings": { "max": 6 },
        "fraction_csp3": { "min": 0.2 }
      },
      {
        "type": "smarts_custom_filter", # reject molecules matching any of these SMARTS
        "patterns": ["[N+](=O)[O-]", "C(=O)Cl"]
      },
      {
        "type": "smarts_catalog_filter", # reject by a named alert catalog
        "catalog": "PAINS" # PAINS | PAINS_A | PAINS_B | PAINS_C | BRENK | CHEMBL | CHEMBL_BMS | CHEMBL_Dundee | CHEMBL_Glaxo | CHEMBL_Inpharmatica | CHEMBL_LINT | CHEMBL_MLSMR | CHEMBL_SureChEMBL | NIH
      },
      {
        "type": "smiles_regex_filter", # reject molecules whose SMILES matches any of these regexes
        "patterns": ["P", "S(=O)(=O)Cl"]
      }
    ]
  }
num_molecules =  100 # how many molecules to generate (10 to 1,000,000)

run_dir = client.small_molecule.design.run(target=target, molecule_filters=molecule_filters, num_molecules=num_molecules, name="test")

