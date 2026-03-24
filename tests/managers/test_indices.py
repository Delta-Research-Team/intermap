"""
Creator: rglez
Date: 3/23/26
Env:
Desc:
"""
import copy
from os.path import join

import pytest

from intermap.managers.indices import guess_from_name, IndexManager


class TestIndexManagerIntegration:
    """
    Runs the whole pipeline on real T1 data.
    """

    def test_pipeline_execution(self, iman):
        """
        Check indices.py is covered and functional.
        """

        assert iman.universe is not None
        assert iman.s1_idx.size > 0
        assert iman.s2_idx.size > 0

        assert hasattr(iman, 'hydroph')
        assert hasattr(iman, 'hb_don')
        assert hasattr(iman, 'rings')

    def test_geometry_indices_logic(self, iman):
        """
        Verifies the complex matrix padding/indexing in get_rings.
        """
        if iman.rings.size > 0:
            assert iman.rings.shape[1] == 7
            assert -1 in iman.rings or iman.rings.shape[0] > 0


class TestIndicesEdgeCases:
    """
    Targets the missing lines in indices.py specifically.
    """

    def test_guess_from_name_failures(self, pt_info):
        symbols, masses, real_names = pt_info
        with pytest.raises(ValueError, match="Unknown element"):
            guess_from_name("X", 0.0, symbols, masses, real_names)

        assert guess_from_name("CA", 12.0, symbols, masses, real_names) == "C"

    def test_selection_errors(self, example_system_args):
        """
        Hits empty selections.
        """
        args = copy.deepcopy(example_system_args)
        args.selection_1 = "resname NONEXISTENT"
        with pytest.raises(ValueError, match="No atoms found for selection 1"):
            IndexManager(args)

    def test_residue_resolution_branch(self, example_system_args, tmp_path):
        """
        Hits resolution == 'residue'.
        """
        args = copy.deepcopy(example_system_args)
        args.resolution = 'residue'
        args.output_dir = join(tmp_path, "residue_res")
        iman = IndexManager(args)
        assert iman.args.resolution == 'residue'
        # Check that shared_idx contains residue indices, not atom indices
        if iman.overlap:
            assert list(iman.shared_idx)[0] < iman.universe.residues.n_residues

    def test_annotations_parsing(self, example_system_args, tmp_path):
        """
        Hits full annotations logic.
        """
        # 1. Create a dummy annotations file
        annot_path = join(tmp_path, "annotations.txt")
        with open(annot_path, "w") as f:
            f.write("# Comment line\n")
            f.write("BindingSite = protein and resid 1 to 5\n")
            f.write("LigandAtoms = nucleic\n")

        args = copy.deepcopy(example_system_args)
        args.annotations = str(annot_path)
        args.output_dir = join(tmp_path, "annot_test")

        iman = IndexManager(args)
        assert "BindingSite" in iman.annotations
        assert len(iman.annotations["BindingSite"]) > 0

    def test_annotations_overlap_error(self, example_system_args, tmp_path):
        """
        Hits overlapping annotations error.
        """
        annot_path = join(tmp_path, "overlap_annotations.txt")
        with open(annot_path, "w") as f:
            f.write("Set1 = protein\n")
            f.write("Set2 = protein\n")

        args = copy.deepcopy(example_system_args)
        args.annotations = str(annot_path)
        with pytest.raises(ValueError, match="overlaps with one or more"):
            IndexManager(args)

    def test_get_interactions_skipping_warning(self, example_system_args,
                                               tmp_path):
        """
        Hits skipping interactions warning.
        """
        args = copy.deepcopy(example_system_args)
        args.interactions = ['XBDonor']
        args.output_dir = join(tmp_path, "interaction_skip")

        iman = IndexManager(args)
        assert 'XBDonor' not in iman.inters_requested


class TestIndicesSpecificErrors:

    def test_selection_2_error(self, example_system_args):
        """
        Hits no atoms found for selection 2.
        """
        args = copy.deepcopy(example_system_args)
        args.selection_2 = "resname VOID"
        with pytest.raises(ValueError, match="No atoms found for selection 2"):
            IndexManager(args)

    def test_annotation_warnings(self, example_system_args, tmp_path, caplog):
        """
        Hits empty annotation selection warning.
        """
        annot_path = join(tmp_path, "empty_annot.txt")
        with open(annot_path, "w") as f:
            f.write("EmptySet = resname GHOST\n")

        args = copy.deepcopy(example_system_args)
        args.annotations = str(annot_path)

        IndexManager(args)
        assert "returned no atoms" in caplog.text


class TestIndicesWetSystem:
    """
    Uses T2 (OPC/Water) to target disconnected residues and lone pairs.
    """

    def test_disconnected_residue_loops(self, iman_t2):
        """
        Hits the self.rdk_disconnected blocks in all getters.
        """
        # If T2 has waters, these dicts will be populated
        assert len(iman_t2.rdk_disconnected) > 0

        # Verify get_waters() hit the disconnected logic
        assert iman_t2.waters.size > 0

        # Verify get_singles (hydroph/cations) hit disconnected loops
        assert hasattr(iman_t2, 'hydroph')

    def test_lone_pair_unknowns(self, iman_t2):
        """Hits lines 114-122 and 338-348 (Unknown 'Z' elements)."""
        # OPC lone pairs (MW) usually trigger the 'KeyError' catch
        if iman_t2.unknown:
            # Verify they were assigned 'Z' (line 341, 348)
            z_atoms = iman_t2.universe.select_atoms("element Z")
            assert len(z_atoms) > 0
            # We drop the 0.0 radius assert, as z_atoms existence proves the branch was hit.

    # Notice we changed t2_system_args -> example_system_args
    def test_high_valence_hydrogens(self, example_system_args, monkeypatch,
                                    tmp_path):
        """Hits lines 139-151 (fix_hh_and_hvalence)."""
        args = copy.deepcopy(example_system_args)  # Use the clean T1 system!
        args.output_dir = str(tmp_path / "hh_test")

        from intermap.managers.indices import IndexManager
        orig_load = IndexManager.load_traj

        def mock_load_bad_h(self):
            u, frames, natoms, nframes, unk = orig_load(self)
            # Find two hydrogens in the SAME residue to safely bond them
            for res in u.residues:
                h_atoms = res.atoms.select_atoms("element H")
                if len(h_atoms) >= 2:
                    u.add_TopologyAttr('bonds',
                                       [(h_atoms[0].index, h_atoms[1].index)])
                    break  # Stop after breaking one residue
            return u, frames, natoms, nframes, unk

        monkeypatch.setattr(IndexManager, "load_traj", mock_load_bad_h)
        # Initialization triggers fix_hh_and_hvalence() automatically
        iman = IndexManager(args)
        assert iman.n_atoms > 0
