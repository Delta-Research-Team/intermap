# Created by rglez at 1/11/25
from os.path import join

import MDAnalysis as mda
import prolif as plf


def run_prolif_test(out_dir):
    # Load universe
    u = mda.Universe(plf.datafiles.TOP, plf.datafiles.TRAJ)

    # Make selections
    sel1 = "protein"
    sel2 = "protein"
    small_protein_selection = u.select_atoms(sel1)
    large_protein_selection = u.select_atoms(sel2,
                                             peptide=small_protein_selection)

    # Run on default interactions
    fp = plf.Fingerprint(count=True)
    fp.run(u.trajectory[::10], small_protein_selection, large_protein_selection)

    # Save
    pickle_name = join(out_dir, "prolif_default.pkl")
    fp.to_pickle(pickle_name)
    return pickle_name, fp.to_dataframe(), fp


# =============================================================================
# Generate the prolif dataframe
# =============================================================================
import time
stamp = time.time()
out_dir = '/home/rglez/RoyHub/intermap/tests/prolif_outputs'
prolif_pkl, df, fp = run_prolif_test(out_dir)
df.columns = ["-".join(a) for a in df.columns.to_flat_index()]
df2 = df.T
plif_time = time.time() - stamp
print("Prolif time:", plif_time)
# =============================================================================
#
# =============================================================================
import rgpack.generals as gnl
imap_pickle = '/home/rglez/RoyHub/intermap/tests/imap_outputs/8oxoGA2_1_InterMap.pickle'
imap_dict = gnl.unpickle_from_file(imap_pickle)
