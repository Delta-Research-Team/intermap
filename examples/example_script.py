# Created by rglez at 12/4/25
from os.path import join

import bitarray.util as bu
import pandas as pd

from intermap import runner

# =============================================================================
# Run InterMap
# =============================================================================

# ----| Define the root directory for InterMap (adjust this path as needed)
my_rootdir = '/home/rglez/RoyHub/intermap/'

# ----| Execute InterMap with the example configuration file
cfg_path = join(my_rootdir, 'examples', 'T1', 'prot-dna_InterMap.cfg')
all_inters = runner.execute(cfg_path=cfg_path)

# ----| Uncompress the bitarrays (only if water bridges were requested)
first_element = all_inters[list(all_inters.keys())[0]]
is_bytes = isinstance(first_element, bytes)
all_inters = {k: bu.sc_decode(v)
              for k, v in all_inters.items()} if is_bytes else all_inters

# =============================================================================
# Create a DataFrame to store the interaction details and their prevalence
# =============================================================================
data = []
for key, values in all_inters.items():
    # Extract the detailed  information from the key
    s1, anot1, s2, annot2, wat, inter_name = key
    s1_resname, s1_resnum, s1_resindex, s1_atomname, s1_atomindex = s1.split(
        '_')
    s2_resname, s2_resnum, s2_resindex, s2_atomname, s2_atomindex = s2.split(
        '_')

    # Compute the prevalence of the interaction
    prevalence = values.count(True) / len(values) * 100

    # Append the interaction details and prevalence to the DataFrame
    df = data.append({
        's1_resname': s1_resname,
        's1_resnum': s1_resnum,
        's1_resindex': s1_resindex,
        's1_atomname': s1_atomname,
        's1_atomindex': s1_atomindex,
        'anot1': anot1,
        's2_resname': s2_resname,
        's2_resnum': s2_resnum,
        's2_resindex': s2_resindex,
        's2_atomname': s2_atomname,
        's2_atomindex': s2_atomindex,
        'anot2': annot2,
        'wat': wat,
        'inter_name': inter_name,
        'prevalence': prevalence
    })

df = pd.DataFrame(data)

# =============================================================================
# From here, you can perform various analyses on the DataFrame `df`, such as:
# - Filtering interactions based on prevalence thresholds
# - Grouping interactions by type or residue
# - Visualizing the distribution of interaction types
# =============================================================================
