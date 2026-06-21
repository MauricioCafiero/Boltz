import os, re, sys, base64
from boltz_api import Boltz
from functions_boltz import get_job_name, make_contraints_template, parse_input

boltz_key = os.getenv('BOLTZ_API_KEY')
client = Boltz(api_key=boltz_key)
print('Read Boltz API key')

job_name=get_job_name()

protein_sequence, ligand_smiles, contact_residues, max_distance,close_residues,\
res_distance, template_type, template_pdbid, template_file = parse_input(job_name)

constraints, template_structure = make_contraints_template(template_type, template_pdbid, template_file, contact_residues, 
                                                        max_distance, close_residues, res_distance)

prediction_input = {
    "entities": [
        {"type": "protein", "value": protein_sequence, "chain_ids": ["A"]},
        {"type": "ligand_smiles", "value": ligand_smiles, "chain_ids": ["B"]},
    ],

    "binding": {"type": "ligand_protein_binding", "binder_chain_id": "B"},

    "num_samples": 3,
}

# Add optional constraints only if they exist
if constraints:
    prediction_input["constraints"] = constraints

# Add optional templates only if they exist
if template_structure:
    prediction_input["templates"] = [
        {
            "template_structure": template_structure,
            "template_chains": [
                { "input_chain_id": "A", "template_chain_id": "A" }
            ],
            "force_threshold_angstroms": 5.0
        }
    ]

'''perform the Boltz prediction via the API using the protein sequence and ligand smiles,
then download the results to a run directory. The run directory will be named 
after the job name'''

# One call: submit, wait, and download the result to a run directory.
submitted_job_name = job_name.replace('.txt', '')
run_dir = client.predictions.structure_and_binding.run(
    model="boltz-2.1", input=prediction_input, name=submitted_job_name
)

print(f'Boltz job completed.')
print('----------------------------------------------------')

