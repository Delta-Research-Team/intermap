# Created by rglez at 5/10/25

from bitarray import bitarray as ba, util as bu

num_frames = 10000
num_inters = 1000
buffer = 100


def iter_decode(bits):
    """
    Decode the bits
    """
    for x in bu.sc_decode(bits):
        yield x


dict_comp = {0: bu.sc_encode(ba('01110')),
             1: bu.sc_encode(ba('11001')),
             2: bu.sc_encode(ba('11111')),
             3: bu.sc_encode(ba('00010'))}

dict_iters = {}
for x in dict_comp:
    dict_iters[x] = iter_decode(dict_comp[x])
traspose = zip(*dict_iters.values())
