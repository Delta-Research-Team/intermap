# Created by rglez at 1/19/25
from rgpack import generals as gnl
from bitarray.util import sc_decode

plif_pickle = '/home/rglez/RoyHub/intermap/tests/outputs/lig_prot_ifp.pkl'
imap_pickle = '/home/rglez/RoyHub/intermap/tests/outputs/imap_lig-prot/lig-prot_InterMap.pickle'

plif_data = gnl.unpickle_from_file(plif_pickle).to_dataframe()
imap_data = gnl.unpickle_from_file(imap_pickle)

selected = {}
for (s1, s2, itype) in imap_data:
    if s1.startswith('LIG_1') and s2.startswith('TYR_359'):
        time = imap_data[(s1, s2, itype)]['time']
        imap_data[(s1, s2, itype)]['time'] = sc_decode(time)
        selected[(s1, s2, itype)] = imap_data[(s1, s2, itype)]
