import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import requests
import os
import json
from tqdm import tqdm
from Bio.PDB import PDBParser

data = pd.read_csv('data/peptide_dd_affinity_sta.csv')
data_filt = data[['Protein.Group', 'is.changed', 'from', 'to']]
all_unique_proteins = data_filt['Protein.Group'].unique()
print('Analyzing {} proteins'.format(len(all_unique_proteins)))

final_df = []

for int_prot in tqdm(all_unique_proteins):

    # download data if it doesn't exist
    dest_json = f'data/protein_data/{int_prot}.json'
    if not os.path.exists(dest_json):
        url = f'https://www.uniprot.org/uniprotkb/{int_prot}.json'
        r = requests.get(url)
        with open(dest_json, 'w') as f:
            f.write(r.text)

    dest_pdb = f'data/protein_data/{int_prot}.pdb'
    if not os.path.exists(dest_pdb):
        url = f'https://alphafold.ebi.ac.uk/files/AF-{int_prot}-F1-model_v4.pdb'
        r = requests.get(url)
        with open(dest_pdb, 'w') as f:
            f.write(r.text)

    # get binding sites
    uniprot_info = json.load(open(dest_json))
    features = uniprot_info['features']
    try:
        binding_sites = [x for x in features if x['type'] == 'Binding site']
        binding_locations = [(x['location']['start']['value'], x['location']['end']['value'], x['ligand']['name']) for x in binding_sites]
    except:
        continue

    if len(binding_locations) == 0:
        continue

    # get AlphaFold predicted structure
    parser = PDBParser()
    structure = parser.get_structure('protein', dest_pdb)
    try:
        model = structure[0]
    except:
        print('Encounter error with protein {}'.format(int_prot))
        continue
    plddts = []
    for chain in model:
        for residue in chain:
            if 'CA' in residue:
                plddt = residue['CA'].get_bfactor()
                plddts.append(plddt)
    plddts = np.array(plddts)
    mean_plddt = np.mean(plddts)

    protein_df = data_filt[data_filt['Protein.Group'] == int_prot].copy()
    protein_df['mean_plddt'] = mean_plddt
    out_df = []

    for start, end, ligand in binding_locations:
        toappend_df = protein_df.copy()
        binding_id = f'{ligand}_{start}-{end}'
        binding_mid = round((start + end) / 2)

        dist_CA_mid = []
        dist_CA_start = []
        dist_CA_end = []
        dist_min_atom_start = []
        dist_min_atom_end = []

        for i, row in protein_df.iterrows():
            try:
                pep_start = int(row['from'])
                pep_end = int(row['to'])
                pep_mid = round((pep_start + pep_end) / 2)

                # Binding site CA
                res_binding = model['A'][binding_mid]
                coord_binding = res_binding['CA'].get_coord()

                # Midpoint CA distance
                res_mid = model['A'][pep_mid]
                dist_CA_mid.append(np.linalg.norm(coord_binding - res_mid['CA'].get_coord()))

                # Start CA distance
                res_start = model['A'][pep_start]
                dist_CA_start.append(np.linalg.norm(coord_binding - res_start['CA'].get_coord()))

                # End CA distance
                res_end = model['A'][pep_end]
                dist_CA_end.append(np.linalg.norm(coord_binding - res_end['CA'].get_coord()))

                # Min atom distance at start
                atoms_start = [atom for atom in res_start]
                dists_start = [np.linalg.norm(coord_binding - atom.get_coord()) for atom in atoms_start]
                dist_min_atom_start.append(min(dists_start))

                # Min atom distance at end
                atoms_end = [atom for atom in res_end]
                dists_end = [np.linalg.norm(coord_binding - atom.get_coord()) for atom in atoms_end]
                dist_min_atom_end.append(min(dists_end))

            except KeyError:
                dist_CA_mid.append(np.nan)
                dist_CA_start.append(np.nan)
                dist_CA_end.append(np.nan)
                dist_min_atom_start.append(np.nan)
                dist_min_atom_end.append(np.nan)

        toappend_df['CA_distance_mid'] = dist_CA_mid
        toappend_df['CA_distance_start'] = dist_CA_start
        toappend_df['CA_distance_end'] = dist_CA_end
        toappend_df['min_atom_distance_start'] = dist_min_atom_start
        toappend_df['min_atom_distance_end'] = dist_min_atom_end
        toappend_df['binding_site_id'] = binding_id

        out_df.append(toappend_df)

    out_df = pd.concat(out_df)
    final_df.append(out_df)

final_df = pd.concat(final_df)
final_df.to_csv('distances.csv', index=False)