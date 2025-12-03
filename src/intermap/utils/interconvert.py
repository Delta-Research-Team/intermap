# Created by rglez at 7/14/25
import sys

from bitarray import util as bu
from rgpack import generals as gnl


def convert_to_csv(pickle, output_csv):
    """
    Convert a dictionary of interactions to a CSV file.

    Args:
        pickle (str): Path to the input pickle containing the interaction data.
        output_csv (str): Path to the output CSV file.
    """
    # Check if the input file is a .pickle file
    if not pickle.endswith('.pickle'):
        print(f'Error: The input file must be a .pickle file, got {pickle}')
        sys.exit(1)

    # Check if the output file is a .csv file
    if not output_csv.endswith('.csv'):
        print(f'Error: The output file must be a .csv file, got {output_csv}')
        sys.exit(1)

    # Open the output CSV file
    bit_dict = gnl.unpickle_from_file(pickle)

    with open(output_csv, 'w') as f:
        f.write('s1,note1,s2,note2,s3,inter_name,prevalence,time\n')
        for key, value in bit_dict.items():
            if isinstance(value, bytes):
                value = bu.sc_decode(value)

            prevalence = f'{value.count(True) / len(value) * 100:3.4f}'
            time = value.to01()
            s1, note1, s2, note2, s3, inter_name = key
            line = f'{s1},{note1},{s2},{note2},{s3},{inter_name},{prevalence},{time}\n'
            f.write(line)


def main():
    """
    Main function to run the conversion from pickle to CSV.
    """
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print('Error: Invalid number of arguments. Expected 2 arguments.')
        print(f'Usage: interconvert <input_pickle_file> <output_csv_file>')
        sys.exit(1)

    # Get the input pickle file and output CSV file from command line arguments
    pickle = sys.argv[1]
    output_csv = sys.argv[2]
    convert_to_csv(pickle, output_csv)
    print(f'Converted {pickle} to {output_csv}')

# pickle = '/media/rglez/Roy5T/RoyData/intermap/1000Frames/NSP13/outs_imap/residue-12-50_InterMap.pickle'
# output_csv = '/media/rglez/Roy5T/RoyData/intermap/1000Frames/NSP13/outs_imap/residue-12-50_InterMap.csv'
