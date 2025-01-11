# Created by rglez at 1/11/25
from os.path import join

import MDAnalysis as mda
import prolif as plf


def run_prolif_test(out_dir):
    # Load universe
    u = mda.Universe(plf.datafiles.TOP, plf.datafiles.TRAJ)

    # Make selections
    sel1 = "resid 119:152"
    sel2 = "protein and not group peptide"
    small_protein_selection = u.select_atoms(sel1)
    large_protein_selection = u.select_atoms(sel2,
                                             peptide=small_protein_selection)

    # Run on default interactions
    fp = plf.Fingerprint([
        "HBDonor",
        "HBAcceptor",
        "PiStacking",
        "PiCation",
        "CationPi",
        "Anionic",
        "Cationic",
    ], count=True)
    fp.run(u.trajectory, small_protein_selection, large_protein_selection)

    # Save
    pickle_name = join(out_dir, "prolif_default.pkl")
    fp.to_pickle(pickle_name)
    return pickle_name, fp.to_dataframe()


# =============================================================================
# Generate the prolif dataframe
# =============================================================================
out_dir = '/home/rglez/RoyHub/intermap/tests/prolif_outputs'
prolif_pkl, df = run_prolif_test(out_dir)
df.columns = ["-".join(a) for a in df.columns.to_flat_index()]
df = df.T

# =============================================================================
#
# =============================================================================
