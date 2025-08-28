# Created by gonzalezroy at 4/2/25
import matplotlib as mpl

cmap1 = 'tab20c'

alpha1 = 0.75
alpha2 = 0.50
alpha3 = 0.25

c1 = 'tab:red'
c2 = 'tab:gray'
c3 = 'tab:blue'

lw1 = 0.5
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

    mpl.rc('axes', labelsize=20)
    mpl.rc('lines', linewidth=8, color='k')
    mpl.rc('font', family='monospace', size=20, monospace='Ubuntu Mono')
    mpl.rc('grid', alpha=0.5, color='gray', linewidth=2, linestyle='--')


def reset_matplotlib():
    """
    Reset matplotlib parameters to default.
    """
    mpl.rcParams.update(mpl.rcParamsDefault)
