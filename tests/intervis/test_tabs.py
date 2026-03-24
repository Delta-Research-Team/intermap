import os

import numpy as np
import pandas as pd
import pytest

from intermap.intervis.app.icsv import CSVFilter
from intermap.intervis.tabs.basetab import BaseTab
from intermap.intervis.tabs.heatmap import HeatMap
from intermap.intervis.tabs.Tab_1 import HeatmapPlot as Tab1HeatmapPlot
from intermap.intervis.tabs.Tab_2 import PrevalencePlot
from intermap.intervis.tabs.Tab_3 import LifetimePlot
from intermap.intervis.tabs.Tab_4 import TimeSeriesPlot
from intermap.intervis.tabs.Tab_5 import InterNetwork, normalize_to_range

# Paths to the real T1 example outputs
EXAMPLES_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', 'examples', 'T1',
                 'outputs'))
REAL_PICKLE = os.path.join(EXAMPLES_DIR, 'prot-dna_InterMap.pickle')
REAL_CFG = os.path.join(EXAMPLES_DIR, 'prot-dna_InterMap.cfg')


# =============================================================================
# FIXTURES
# =============================================================================

@pytest.fixture
def real_df(example_system_args, monkeypatch):
    """Loads the REAL data from the T1 pickle file using the actual CSVFilter."""
    import MDAnalysis as mda

    # Patch MDAnalysis Universe so it finds the topology file from the fixture
    # instead of crashing on relative paths in the .cfg
    orig_univ = mda.Universe

    def mock_univ(topo, *args, **kwargs):
        return orig_univ(example_system_args.topology, *args, **kwargs)

    monkeypatch.setattr(mda, "Universe", mock_univ)

    # Parse the real data
    cf = CSVFilter(REAL_PICKLE, REAL_CFG)

    # Take a subset to keep the test blazing fast, but it's 100% real data
    df = cf.master.head(50).copy()

    # Inject a few edge cases missing from T1 to guarantee 100% branch coverage
    # (T1 has no WaterBridges and timeseries are strictly numpy arrays)
    edge_cases = pd.DataFrame({
        'sel1': ['WAT_1_0_O_0', 'ALA_2_1_N_1', 'ALA_3_2_N_2'],
        'sel2': ['TRP_4_3_O_3', 'CYS_5_4_C_4', 'CYS_6_5_C_5'],
        'interaction_name': ['WaterBridge', 'HBDonor', 'PiStacking'],
        'prevalence': [90.0, 50.0, 10.0],
        'note1': ['wat', '', ''],
        'note2': ['', '', ''],
        'timeseries': [
            np.array([1, 1, 0]),  # Normal array
            (1, 1, 0),  # Tuple format (hits Tab_3 branch)
            '101'  # String format (hits Tab_3 branch)
        ],
        'idx1': [0, 1, 2],
        'idx2': [3, 4, 5],
        'water': ['', '', '']
    })

    df = pd.concat([df, edge_cases], ignore_index=True)
    return df


@pytest.fixture
def empty_df(real_df):
    """Provides an empty DataFrame that retains the exact real column schema."""
    return real_df.iloc[0:0].copy()


# =============================================================================
# TESTS
# =============================================================================

def test_basetab(real_df):
    """BaseTab should raise NotImplementedError if process_data isn't overridden."""
    with pytest.raises(NotImplementedError):
        BaseTab(real_df)


def test_heatmap_and_tab1(real_df, empty_df):
    """Tests both HeatMap implementations (heatmap.py and Tab_1.py)."""
    # 1. Test heatmap.py
    hm = HeatMap(real_df, width=800, height=600, show_prevalence=True)
    fig = hm.create_heatmap_plot('X-Axis', 'Y-Axis')
    assert fig is not None
    assert len(fig.data) > 0
    assert HeatMap(empty_df).create_heatmap_plot('X', 'Y') is None

    # 2. Test Tab_1.py (Legacy/Alternative HeatmapPlot)
    hm1 = Tab1HeatmapPlot(real_df, show_prevalence=True)
    fig1 = hm1.create_heatmap_plot('X-Axis', 'Y-Axis')
    assert fig1 is not None
    assert Tab1HeatmapPlot(empty_df).create_heatmap_plot('X', 'Y') is None


def test_tab2_prevalence(real_df, empty_df):
    """Tests the bar charts for Prevalence."""
    p1 = PrevalencePlot(real_df, plot_type='sel1')
    fig1 = p1.create_prevalence_plot('X-Axis', 'Y-Axis')
    assert fig1 is not None

    p2 = PrevalencePlot(real_df, plot_type='sel2')
    fig2 = p2.create_prevalence_plot('X-Axis', 'Y-Axis')
    assert fig2 is not None

    assert PrevalencePlot(empty_df).create_prevalence_plot('X', 'Y') is None


def test_tab3_lifetime(real_df, empty_df):
    """Tests the box plots for Lifetimes."""
    lp = LifetimePlot(real_df)
    fig = lp.create_lifetime_plot('X-Axis', 'Y-Axis')
    assert fig is not None

    assert LifetimePlot(empty_df).create_lifetime_plot('X', 'Y') is None


def test_tab4_timeseries(real_df, empty_df):
    """Tests the multi-trace subplots for interactions over time."""
    ts = TimeSeriesPlot(real_df)
    fig = ts.create_time_series_plot('X-Axis', 'Y-Axis')
    assert fig is not None

    assert TimeSeriesPlot(empty_df).create_time_series_plot('X', 'Y') is None


def test_tab5_network(real_df, empty_df):
    """Tests the PyVis network graph generation."""
    assert normalize_to_range([]) == []
    assert normalize_to_range([10, 20, 30], target_min=0, target_max=1) == [
        0.0, 0.5, 1.0]

    net = InterNetwork(real_df)
    G = net.get_graph()
    assert len(G.nodes) > 0
    assert len(G.edges) > 0

    pynet = net.create_network_plot()
    assert pynet is not None

    # Check that our injected WAT residue was processed correctly
    assert any("WAT" in str(node['id']) for node in pynet.nodes)

    assert InterNetwork(empty_df).create_network_plot() is None
