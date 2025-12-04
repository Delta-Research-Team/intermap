# Created by rglez at 12/4/25

import bitarray.util as bu

from intermap import runner

# ----| 1. Execute InterMap with the given configuration file
cfg_path = './prot-dna_InterMap.cfg'
bit_dict = runner.execute(cfg_path=cfg_path)
first_element = bit_dict[list(bit_dict.keys())[0]]

# ----| 2. Decode the bitarrays if they are in bytes format
if isinstance(first_element, bytes):
    bit_dict = {k: bu.sc_decode(v) for k, v in bit_dict.items()}
else:
    bit_dict = bit_dict

# ----| 3. Continue with further processing of bit_dict as needed
print(f'Processed {len(bit_dict)} interactions from configuration: {cfg_path}')
