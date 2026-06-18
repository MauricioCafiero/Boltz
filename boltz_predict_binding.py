import os, re, sys
from boltz_api import Boltz
from sequence_bank import sequence_bank, smiles_pattern, smiles_node

boltz_key = os.getenv('BOLTZ_API_KEY')
print('got key')
client = Boltz(api_key=boltz_key)
print('used key')

'''
Jobname should be a text file with the job name in it, e.g. "my_job_name.txt"
structure should be: 
protein: [sequence in one letter format]
ligand_smiles: [SMILES string]

currently only supports one protein chain and one ligand, 
'''

job_name = sys.argv[1] if len(sys.argv) > 1 else None
if job_name is None:
    job_name = input('enter job name: ')

if not job_name.endswith('.txt'):
    job_name += '.txt'

print('----------------------------------------------------')
print(f'Boltz job name: {job_name.replace(".txt", "")}')

with open(job_name, 'r') as f:
    lines = f.readlines()
    protein_sequence = lines[0].strip().split(':')[1].strip()
    ligand_smiles = lines[1].strip().split(':')[1].strip()

#------------------------------------------------------------
if protein_sequence.split(',')[0] == 'bank' and protein_sequence.split(',')[1] in sequence_bank.keys():
    protein_sequence = sequence_bank[protein_sequence.split(',')[1]]
    print(f'Protein sequence found in sequence bank')

#check if ligand_smiles is a valid SMILES string using regex using imported smiles pattern from sequence_bank.py
if not re.match(smiles_pattern, ligand_smiles):
    print(f'getting SMILES string from PubChem for {ligand_smiles}')
    ligand_smiles = smiles_node(ligand_smiles)
    if ligand_smiles == 'unknown':
        print(f'Could not find SMILES string for {ligand_smiles} in PubChem')
        sys.exit(1)

#------------------------------------------------------------

prediction_input = {
    "entities": [
        {"type": "protein", "value": protein_sequence, "chain_ids": ["A"]},
        {"type": "ligand_smiles", "value": ligand_smiles, "chain_ids": ["B"]},
    ],
    "binding": {"type": "ligand_protein_binding", "binder_chain_id": "B"},
    "num_samples": 3,
}  # see Input format for constraints, bonds, templates, model_options, …

# One call: submit, wait, and download the result to a run directory.
submitted_job_name = job_name.replace('.txt', '')
run_dir = client.predictions.structure_and_binding.run(
    model="boltz-2.1", input=prediction_input, name=submitted_job_name
)

print(f'Boltz job completed.')
print('----------------------------------------------------')

