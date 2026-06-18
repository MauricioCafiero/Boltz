import json
import sys

#look for name as a command line argument, if not found, ask for input
name = sys.argv[1] if len(sys.argv) > 1 else None
if name is None:
    name = input('enter job name: ')

print('----------------------------------------------------')
print(f'Boltz results for job: {name}')

path = f'boltz-experiments/{name}/outputs/files/prediction/'

metrics = json.load(open(path + 'metrics.json', 'r'))

# score is in 'binding_metrics' -> 'optimization_score'
score = metrics['binding_metrics']['optimization_score']
binding_confidence = metrics['binding_metrics']['binding_confidence']

print(f'optimization score: {score}')
print(f'binding confidence: {binding_confidence}')

for idx, test_run in enumerate(metrics['all_sample_results']):
    struct_confidence = test_run['metrics']['structure_confidence']
    print(f'structure confidence (run {idx}): {struct_confidence}')

total_runs = idx + 1
print(f'total runs: {total_runs}')

for idx in range(total_runs):
    filename = f'sample_{idx}_predicted_structure.cif'
    #copy file to current directory
    with open(path + filename, 'r') as f:
        new_filename = f'{name}_{idx}_predicted_structure.cif'
        with open(new_filename, 'w') as f_out:
            f_out.write(f.read())

print(f'copied {total_runs} predicted structures to current directory')
print('----------------------------------------------------')
