# Created by gonzalezroy at 11/14/24

import numpy as np


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


def get_coordinates(u, chunk, sel_idx):
    xyz_chunk = np.empty((chunk.size, sel_idx.size, 3), dtype=np.float32)
    for i, frame in enumerate(chunk):
        xyz_chunk[i] = u.trajectory[frame].positions[sel_idx].copy()
    return xyz_chunk
