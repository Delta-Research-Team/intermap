# Created by fajardo at 9/1/25

import csv
import pickle
import os

import numpy as np
from bitarray.util import sc_decode


# =========================================================================
# Creating the interaction_results.csv file
# =========================================================================

pickle_path = "/home/fajardo/03_Fajardo_Hub/02_InterMap/intermap/example/imap_output/config_files/outputs/prolif_tutorial/lig-prot/"

pickle_files = [f for f in os.listdir(pickle_path) if f.endswith('.pickle')]
if not pickle_files:
    raise FileNotFoundError("No .pickle file was found in the specified directory :( "
                            "Please verify whether the simulation concluded successfully.")

pickle_file = os.path.join(pickle_path, pickle_files[0])


print(f"Loading data from: {pickle_file}...")
with open(pickle_file, 'rb') as f:
    all_data = pickle.load(f)


filename = f"interaction_results.csv"

all_decoded_rows = []

for key, value in all_data.items():
    ligand, protein, interaction = key
    time_data = value['time']

    decoded_data = sc_decode(time_data)

    row_data = {
        'ligand': ligand,
        'protein': protein,
        'interaction': interaction,
        'states': [1 if state else 0 for state in decoded_data]
    }
    all_decoded_rows.append(row_data)


print(f"Writing data in {filename}...")
with open(filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)

    num_frames = len(all_decoded_rows[0]['states'])
    headers = ['ligand', 'protein', 'interaction'] + [str(i) for i in range(num_frames)]
    writer.writerow(headers)

    for row in all_decoded_rows:
        data_row = [row['ligand'], row['protein'], row['interaction']] + row['states']
        writer.writerow(data_row)


print(f"File CSV succesfully created: {filename} ;)")
