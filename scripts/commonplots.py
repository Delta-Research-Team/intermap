# Created by gonzalezroy at 4/2/25
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# choose number of discrete colors


N = 10  # number of discrete colors
colors = plt.cm.RdYlGn_r(np.linspace(0.05, 0.95, N))
cmap1 = matplotlib.colors.ListedColormap(colors)

alpha1 = 0.75
alpha2 = 0.50
alpha3 = 0.25

c1 = 'tab:red'
c2 = 'tab:orange'
c3 = 'tab:green'

lw1 = 1
fs1 = 20


# =============================================================================
#
# =============================================================================
def generic_matplotlib():
    """
    Some customizations of matplotlib.
    """
    print('Customizing matplotlib')
    mpl.rc('figure', figsize=[12, 8], dpi=300)
    mpl.rc('xtick', direction='in', top=True)
    mpl.rc('xtick.major', top=False)
    mpl.rc('xtick.minor', top=False, visible=True, bottom=False)
    mpl.rc('ytick', direction='in', right=True)
    mpl.rc('ytick.major', right=True, )
    mpl.rc('ytick.minor', right=True, visible=True)

    mpl.rc('lines', linewidth=8, color='k')
    mpl.rc('grid', alpha=0.5, color='gray', linewidth=2, linestyle='--')

    mpl.rc('font', family='STIXGeneral', size=22, monospace='stix')
    mpl.rc('axes', labelsize=22)
    mpl.rcParams['xtick.labelsize'] = 20  # or whatever size you want
    mpl.rcParams['ytick.labelsize'] = 20

    # mpl.rcParams['mathtext.fontset'] = 'stix'


def reset_matplotlib():
    """
    Reset matplotlib parameters to default.
    """
    mpl.rcParams.update(mpl.rcParamsDefault)
