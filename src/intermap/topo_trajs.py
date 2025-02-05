# Created by gonzalezroy at 11/14/24
import os
from collections import defaultdict
from os.path import join

import numpy as np
from numba.typed.typeddict import Dict


def get_resids_indices(trajectory):
    """
    Get indices of residues in the load trajectory

    Args:
        trajectory: trajectory loaded in mdtraj format

    Returns:
        res_ind_numba: numba Dict of each residue's indices
        babel_dict: the equivalence between the original resid numbering and
                    the 0-based used internally
    """
    # Parse the topological information
    df = trajectory.topology.to_dataframe()[0]

    group_by_index = df.groupby(["chainID", "resSeq", "segmentID"]).indices
    group_by_index_noh = {}
    for key in group_by_index:
        values = group_by_index[key]
        noh = values[df.loc[values, "element"] != "H"]
        group_by_index_noh[key] = noh

    babel_dict = {i: x for i, x in enumerate(group_by_index)}

    # Transform to numba-dict
    res_ind_zero = {i: group_by_index[x] for i, x in enumerate(group_by_index)}
    res_ind_noh = {i: group_by_index_noh[x] for i, x in
                   enumerate(group_by_index_noh)}

    res_ind_numba = pydict_to_numbadict(res_ind_zero)
    res_ind_noh_numba = pydict_to_numbadict(res_ind_noh)
    return res_ind_numba, res_ind_noh_numba, babel_dict


def pydict_to_numbadict(py_dict):
    """
    Converts from Python dict to Numba dict

    Args:
        py_dict: Python dict

    Returns:
        numba_dict: Numba dict
    """
    numba_dict = Dict()
    for key in py_dict:
        numba_dict.update({key: py_dict[key]})
    return numba_dict


def get_dha_indices(trajectory, heavies_elements, atoms_to_resids):
    """
    Get 0-based indices of donors, hydrogens, and acceptors in an MDTraj traj

    Args:
        trajectory: MDTraj trajectory object
        heavies_elements: name of elements considered as heavies
        atoms_to_resids: dict mapping atoms indices to residues indices

    Returns:
        donors: indices of donor atoms (N or O bonded to H)
        hydros: indices of hydrogen atoms (H bonded to N or O)
        heavies: indices of heavy atoms (N or O)
    """
    # Get heavies and hydrogen indices
    df, bonds = trajectory.topology.to_dataframe()
    all_hydrogens = set(df[df.element == "H"].index)

    a_raw1 = set(np.where(df.element.isin(heavies_elements, ))[0])

    # Find D-H indices
    h_raw1 = []
    d_raw1 = []
    for values in bonds:
        at1 = int(values[0])
        at2 = int(values[1])
        if (at1 in all_hydrogens) and (at2 in a_raw1):
            h_raw1.append(at1)
            d_raw1.append(at2)
        elif (at2 in all_hydrogens) and (at1 in a_raw1):
            h_raw1.append(at2)
            d_raw1.append(at1)
        else:
            continue

    # Process the indices of donors to a numba dict
    d_raw = defaultdict(list)
    [d_raw[atoms_to_resids[x]].append(x) for x in d_raw1]
    d_raw3 = {x: np.asarray(d_raw[x], dtype=np.int32) for x in d_raw}
    donors = pydict_to_numbadict(d_raw3)

    # Process the indices of hydrogens to a numba dict
    h_raw = defaultdict(list)
    [h_raw[atoms_to_resids[x]].append(x) for x in h_raw1]
    h_raw3 = {x: np.asarray(h_raw[x], dtype=np.int32) for x in h_raw}
    hydros = pydict_to_numbadict(h_raw3)

    # Process the indices of acceptors to a numba dict
    a_raw = defaultdict(list)
    [a_raw[atoms_to_resids[x]].append(x) for x in a_raw1]
    a_raw3 = {x: np.asarray(a_raw[x], dtype=np.int32) for x in a_raw}
    acceptors = pydict_to_numbadict(a_raw3)
    return donors, hydros, acceptors


def mount_data_dir(proj_dir, env_var):
    """
    Mount the data directory to the project directory
    """
    data_dir = os.getenv(env_var)

    # Check if the environment variable is defined
    if not data_dir:
        raise ValueError(f'The environment variable {env_var} is not set. '
                         'Please define it at /etc/profile.')

    # Check if the directory defined by the environment variable exists
    if not os.path.exists(data_dir):
        raise FileNotFoundError(
            f'The directory defined by {env_var} does not exist: {data_dir}')

    # Define the local data directory
    local_data_dir = join(proj_dir, 'data')
    os.unlink(local_data_dir) if os.path.exists(local_data_dir) else None
    os.symlink(data_dir, local_data_dir)
    return local_data_dir


def parse_last_param(last, traj_len):
    """
    Parse the last parameter

    Args:
        last (int): Last frame to process
        traj_len (int): Length of the trajectory

    Returns:
        last (int): Last frame to process
    """
    if last == -1:
        last = traj_len
    elif last > traj_len:
        print(
            f"WARNING: The declared end frame {last} is larger than the trajectory "
            f"length {traj_len}. The end frame will be set to the last frame.")
        last = traj_len
    return last


def split_in_chunks(array, chunk_size):
    """
    Split an array in chunks of the specified size

    Args:
        array (ndarray): Array to split
        chunk_size (int): Size of the chunks

    Returns:
        chunks (generator): Generator with the chunks
    """
    for i in range(0, array.shape[0], chunk_size):
        yield array[i:i + chunk_size]


def get_coordinates(u, chunk, sel_idx, num_atoms):
    xyz_chunk = np.empty((chunk.size, num_atoms, 3), dtype=np.float32)
    for i, frame in enumerate(chunk):
        xyz_chunk[i] = u.trajectory[frame].positions[sel_idx]
    return xyz_chunk
