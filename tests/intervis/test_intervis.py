import os

import MDAnalysis as mda
import numpy as np
import pandas as pd
import pytest

from intermap.intervis.app.icsv import (compress_wb, CSVFilter,
                                        process_heatmap_data,
                                        process_lifetime_data,
                                        process_prevalence_data,
                                        process_prevalence_data2,
                                        process_time_series_data, sortby,
                                        transpose)

# =============================================================================
# FIXTURES
# =============================================================================

EXAMPLES_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', 'examples', 'T1',
                 'outputs'))
REAL_PICKLE = os.path.join(EXAMPLES_DIR, 'prot-dna_InterMap.pickle')
REAL_CFG = os.path.join(EXAMPLES_DIR, 'prot-dna_InterMap.cfg')


@pytest.fixture
def real_csvf(example_system_args, monkeypatch):
    """Loads the real T1 data. Patches Universe to ensure the PDB is found."""
    orig_univ = mda.Universe

    def mock_univ(topo, *args, **kwargs):
        # Force the correct absolute path from your existing fixture
        return orig_univ(example_system_args.topology, *args, **kwargs)

    monkeypatch.setattr(mda, "Universe", mock_univ)
    return CSVFilter(REAL_PICKLE, REAL_CFG)


@pytest.fixture
def master_df(real_csvf):
    return real_csvf.master.copy()


# =============================================================================
# DATA MANIPULATION TESTS
# =============================================================================

def test_transpose(master_df):
    df_t = transpose(master_df.copy())
    assert 'sel1' in df_t.columns
    # Check that columns swapped correctly by verifying original sel2 is now sel1
    assert df_t['sel1'].iloc[0] == master_df['sel2'].iloc[0]


def test_sortby(master_df):
    df1 = sortby(master_df, 'name')
    df2 = sortby(master_df, 'number')
    df3 = sortby(master_df, 'annotation')

    assert not df1.empty and not df2.empty and not df3.empty

    with pytest.raises(ValueError, match="Invalid choice"):
        sortby(master_df, 'invalid')


def test_compress_wb():
    """Injects a WaterBridge to test the compression logic since T1 lacks them."""
    df = pd.DataFrame({
        'sel1': ['A_1_0_N_0', 'A_1_0_N_0'],
        'sel2': ['C_3_2_C_2', 'C_3_2_C_2'],
        'interaction_name': ['WaterBridge', 'WaterBridge'],
        'timeseries': ['100', '010']
    })
    df_comp = compress_wb(df.copy())

    assert len(df_comp) == 1
    assert np.array_equal(df_comp['timeseries'].iloc[0], [1, 1, 0])


# =============================================================================
# CSVFILTER METHOD TESTS
# =============================================================================

def test_csvfilter_filters(real_csvf):
    # Prevalence
    idx, stat = real_csvf.by_prevalence(0.0)
    assert stat == 0 and len(idx) > 0
    idx, stat = real_csvf.by_prevalence(101.0)
    assert stat == -1

    # Interactions
    idx, stat = real_csvf.by_inters('all')
    assert stat == 0
    inter_name = list(real_csvf.inters2df.keys())[0]
    idx, stat = real_csvf.by_inters([inter_name, 'FakeInter'])
    assert stat == 0

    # Notes
    idx, stat = real_csvf.by_notes('all')
    assert stat in (0, -1)  # T1 might not have annotations
    if real_csvf.notes2df:
        note_name = list(real_csvf.notes2df.keys())[0]
        idx, stat = real_csvf.by_notes([note_name])
        assert stat == 0

    # MDA Selections
    idx, stat = real_csvf.by_mda('all')
    assert stat == 0
    idx, stat = real_csvf.by_mda('resname FAKENAME')
    assert stat == -1


# =============================================================================
# PLOT PROCESSORS TESTS
# =============================================================================

def test_plot_processors(master_df):
    df = master_df.copy()

    # process_prevalence_data2 requires these columns to be parsed beforehand
    df['resname1'] = 'ALA'
    df['resnum1'] = 1
    df['resname2'] = 'GLY'
    df['resnum2'] = 2

    # Heatmap Data
    hm = process_heatmap_data(df)
    assert 'pivot_interaction' in hm

    # Prevalence Data
    pd1 = process_prevalence_data(df, 'sel1')
    assert isinstance(pd1, list)

    # Prevalence Data 2 (Multiple Sorts)
    for sort_col in ['resname', 'resnum', 'idx', 'note']:
        res = process_prevalence_data2(df.copy(), 'sel1', sort_by=sort_col)
        assert isinstance(res, list)

    with pytest.raises(ValueError, match="Invalid sort_by"):
        process_prevalence_data2(df, 'sel1', sort_by='invalid')

    # Time Series Data
    ts = process_time_series_data(df)
    assert 'scatter_df' in ts


def test_lifetime_data_formats(master_df):
    """Forces multiple data types to hit all formatting branches in lifetime."""
    df = master_df.head(4).copy()

    # 1. Standard Array/List (from real data) -> yields 2 intervals [0,1] and [2,3]
    df.at[0, 'timeseries'] = [1, 0, 1]

    # 2. Tuple format -> yields 1 interval [0,2]
    df.at[1, 'timeseries'] = (1, 1, 0)

    # 3. String format -> yields 2 intervals [0,2] and [3,4]
    df.at[2, 'timeseries'] = '1101'

    # 4. Unsupported format (triggers print and skip) -> yields 0 intervals
    df.at[3, 'timeseries'] = 123.45

    res = process_lifetime_data(df)
    assert not res.empty

    # Total intervals = 2 + 1 + 2 + 0 = 5
    assert len(res) == 5


def test_icsv_edge_cases(tmp_path, example_system_args, monkeypatch):
    """Hits remaining branches: prevalence 0/100, and empty note filtering."""
    import pickle
    import bitarray.util as bu
    from intermap.intervis.app.icsv import parse_pickle, CSVFilter

    # 1. Test parse_pickle with uncompressed bitarrays (hits line 140)
    p_path = tmp_path / "uncompressed.pickle"
    # Create a 4-frame bitarray: [1, 1, 1, 1] (100% prevalence)
    b_full = bu.zeros(4)
    b_full.setall(1)
    d = {('ALA_1_0_N_0', '', 'GLY_2_1_CA_1', '', '', 'Contact'): b_full}
    with open(p_path, 'wb') as f:
        pickle.dump(d, f)

    df, res = parse_pickle(str(p_path))
    assert df['prevalence'].iloc[0] == 100.0

    # 2. Test CSVFilter notes branch when no notes exist (hits line 355-356)
    # We reuse the REAL_CFG but a dummy pickle with no notes
    c_path = os.path.join(EXAMPLES_DIR, 'prot-dna_InterMap.cfg')

    # Force MDAnalysis to find the topology
    import MDAnalysis as mda
    orig_univ = mda.Universe
    monkeypatch.setattr(mda, "Universe", lambda topo, *args: orig_univ(
        example_system_args.topology))

    cf = CSVFilter(str(p_path), c_path)
    # Target line 356: filter by notes when notes2df is empty
    idx, stat = cf.by_notes(['SomeNote'])
    assert stat == -1
