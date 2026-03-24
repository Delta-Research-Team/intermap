import os

import pandas as pd
import pytest

# Target modules
import intermap.intervis.app.helpers as hp
import intermap.intervis.app.plots as pl


# =============================================================================
# TESTS FOR HELPERS.PY
# =============================================================================

class TestUIHelpers:
    def test_find_topology_file(self, tmp_path):
        """Tests the file extension scanning logic."""
        # Setup dummy directory
        dcd_path = tmp_path / "trajectory.dcd"
        dcd_path.touch()

        # 1. Test Exact Match (same basename)
        pdb_path = tmp_path / "trajectory.pdb"
        pdb_path.touch()

        found, path = hp.find_topology_file(None, str(dcd_path))
        assert found is True
        assert path == str(pdb_path)

        # 2. Test Fallback (different basename, valid extension)
        os.remove(pdb_path)
        gro_path = tmp_path / "other_name.gro"
        gro_path.touch()

        found, path = hp.find_topology_file(None, str(dcd_path))
        assert found is True
        assert path == str(gro_path)

        # 3. Test Failure (No valid extensions)
        os.remove(gro_path)
        found, path = hp.find_topology_file(None, str(dcd_path))
        assert found is False

    def test_get_image_base64(self, tmp_path):
        """Tests image encoding and error handling."""
        img_path = tmp_path / "test.png"
        img_path.write_bytes(b"dummy image data")

        # Test success
        encoded = hp.get_image_base64(str(img_path))
        assert encoded.startswith("data:image/jpeg;base64,")

        # Test failure
        encoded_fail = hp.get_image_base64(str(tmp_path / "nonexistent.png"))
        assert encoded_fail == ""

    def test_calculate_prevalence(self):
        row = pd.Series({'prevalence': 75.5, 'other': 1})
        assert hp.calculate_prevalence(row) == 75.5
        assert hp.calculate_prevalence(pd.Series({'other': 1})) is None

    def test_generate_interaction_choices(self):
        assert hp.generate_interaction_choices(None) == {}

        df = pd.DataFrame(
            {'interaction_name': ['HBDonor', 'HBDonor', 'PiStacking']})
        choices = hp.generate_interaction_choices(df)

        assert 'HBDonor' in choices
        assert choices['HBDonor'] == 'HBDonor (2)'
        assert choices['PiStacking'] == 'PiStacking (1)'

    def test_validate_mda_selection(self):
        # 1. Valid
        valid, msg = hp.validate_mda_selection("resname ALA and backbone")
        assert valid is True

        # 2. Empty
        valid, msg = hp.validate_mda_selection("   ")
        assert valid is False and "empty" in msg

        # 3. Invalid characters
        valid, msg = hp.validate_mda_selection("protein @")
        assert valid is False and "invalid characters" in msg

        # 4. No valid keywords
        valid, msg = hp.validate_mda_selection("random unknown words")
        assert valid is False and "valid MDAnalysis keywords" in msg


# =============================================================================
# TESTS FOR PLOTS.PY (Mocked to avoid rendering)
# =============================================================================

class MockPlotObj:
    """Mocks the behavior of the internal plotting classes."""

    def __init__(self, *args, **kwargs):
        pass

    def create_heatmap_plot(self, *args, **kwargs):
        return "heatmap_fig"

    def create_prevalence_plot(self, *args, **kwargs):
        return "prevalence_fig"

    def create_lifetime_plot(self, *args, **kwargs):
        return "lifetime_fig"

    def create_time_series_plot(self, *args, **kwargs):
        return "timeseries_fig"

    def create_network_plot(self, *args, **kwargs):
        return "network_fig"


class TestPlotRouters:
    @pytest.fixture(autouse=True)
    def patch_plot_classes(self, monkeypatch):
        """Bypasses the actual Plotly/Network graph generation."""
        monkeypatch.setattr(pl, "HeatMap", MockPlotObj)
        monkeypatch.setattr(pl, "PrevalencePlot", MockPlotObj)
        monkeypatch.setattr(pl, "LifetimePlot", MockPlotObj)
        monkeypatch.setattr(pl, "TimeSeriesPlot", MockPlotObj)
        monkeypatch.setattr(pl, "InterNetwork", MockPlotObj)

    def test_empty_df_returns(self):
        """Tests the empty DataFrame escape hatches."""
        empty_df = pd.DataFrame()
        assert pl.create_plot(empty_df, 1, 1, 'x', 'y') is None
        assert pl.create_lifetime_plot(empty_df, 1, 1, 'x', 'y') is None
        assert pl.create_network_plot(empty_df, 1, 1, 'x', 'y') is None

    def test_plot_routers(self):
        """Tests that the router functions successfully call their target classes."""
        df = pd.DataFrame({'dummy': [1]})

        assert pl.create_plot(df, 100, 100, 'x', 'y') == "heatmap_fig"
        assert pl.create_sel1_interactions_plot(df, 100, 100, 'x',
                                                'y') == "prevalence_fig"
        assert pl.create_sel2_interactions_plot(df, 100, 100, 'x',
                                                'y') == "prevalence_fig"
        assert pl.create_lifetime_plot(df, 100, 100, 'x',
                                       'y') == "lifetime_fig"
        assert pl.create_interactions_over_time_plot(df, 100, 100, 'x',
                                                     'y') == "timeseries_fig"

        # Test Network Plot with and without params
        assert pl.create_network_plot(df, 100, 100, 'x', 'y') == "network_fig"
        assert pl.create_network_plot(df, 100, 100, 'x', 'y',
                                      {'gravity': 0}) == "network_fig"


def test_find_topology_exception():
    """Triggers the Exception block in find_topology_file (lines 44-46)."""
    # Passing None triggers a TypeError inside pathlib
    found, path = hp.find_topology_file(None, None)
    assert found is False
    assert path == ""
