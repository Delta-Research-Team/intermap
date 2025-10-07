# Created by gonzalezroy at 5/16/25
import matplotlib.pyplot as plt
from bitarray import bitarray as ba, util as bu

array = ba('00010101110101011111000010101011011111110010111111010101000111')
lifetimes = [x[2] - x[1] for x in bu.intervals(array) if x[0] == 1]
plt.boxplot(lifetimes)
plt.show()
