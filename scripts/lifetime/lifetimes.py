# Created by gonzalezroy at 5/16/25
from bitarray import bitarray as ba
from bitarray import util as bu
import numpy as np
import matplotlib.pyplot as plt
array = ba('00010101110101011111000010101011011111110010111111010101000111')
lifetimes = np.mean([x[2] - x[1] for x in bu.intervals(array) if x[0] == 1])
plt.hist(lifetimes)
plt.show()
