from unittest.mock import MagicMock

import pandas as pd
import plotly.graph_objects as go

from intermap.intervis.app.main import (apply_organization_to_figure,
                                        get_axis_names, process_axis_labels,
                                        simplify_pair_label)


def test_get_axis_names():
    assert get_axis_names(None) == ("sel1", "sel2")

    mock_obj = MagicMock()
    mock_obj.axisx = "Custom X"
    mock_obj.axisy = "Custom Y"
    assert get_axis_names(mock_obj) == ("Custom X", "Custom Y")

    mock_obj.axisx = ""
    mock_obj.axisy = None
    assert get_axis_names(mock_obj) == ("sel1", "sel2")


def test_simplify_pair_label():
    # Valid
    assert simplify_pair_label(
        "ALA_1_N - GLY_2_O (HBDonor)") == "ALA_1 - GLY_2 (HBDonor)"
    assert simplify_pair_label(
        "ALA_1 - GLY_2 (HBDonor)") == "ALA - GLY (HBDonor)"

    # Invalid / Edges
    assert simplify_pair_label("No Separator") == "No Separator"
    assert simplify_pair_label(None) is None
    assert simplify_pair_label(
        "A_1 - B_2 - C_3 (Test)") == "A_1 - B_2 - C_3 (Test)"


def test_process_axis_labels():
    df = pd.DataFrame({
        'sel1': ['ALA_1_N'],
        'sel2': ['GLY_2_O'],
        'pair': ['ALA_1_N - GLY_2_O (Test)'],
        'selection_pair': ['ALA_1_N - GLY_2_O (Test)']
    })

    # Test simplify off
    assert process_axis_labels(df, False) is df

    # Test simplify on
    res = process_axis_labels(df, True)
    assert res['sel1'][0] == "ALA_1"
    assert res['sel2'][0] == "GLY_2"
    assert res['pair'][0] == "ALA_1 - GLY_2 (Test)"
    assert res['selection_pair'][0] == "ALA_1 - GLY_2 (Test)"


def test_apply_organization_to_figure():
    fig = go.Figure(data=[go.Scatter(x=["ALA_2_N", "GLY_1_O__A_B", "TRP"])])
    if fig is None or not getattr(fig, "data", None) or len(
            fig.data) == 0 or not hasattr(fig.data[0], 'x') or fig.data[
        0].x is None:
        return fig

    df = pd.DataFrame(
        {'sel1': ["ALA_2_N", "GLY_1_O__A_B", "TRP"], 'note1': ["C", "A", "B"]})

    # Test resname sort
    fig_resname = apply_organization_to_figure(fig, "resname")
    assert list(fig_resname.layout.xaxis.categoryarray) == ["ALA_2_N",
                                                            "GLY_1_O__A_B",
                                                            "TRP"]

    # Test resnum sort
    fig_resnum = apply_organization_to_figure(fig, "resnum")
    assert list(fig_resnum.layout.xaxis.categoryarray) == ["TRP",
                                                           "GLY_1_O__A_B",
                                                           "ALA_2_N"]

    # Test annotation sort (Requires dataframe)
    fig_annot = apply_organization_to_figure(fig, "annotation", df)
    assert list(fig_annot.layout.xaxis.categoryarray) == ["GLY_1_O__A_B",
                                                          "TRP", "ALA_2_N"]

    # Escape hatches
    assert apply_organization_to_figure(None, "resname") is None
    fig_empty = go.Figure()
    assert apply_organization_to_figure(fig_empty, "resname") is fig_empty
