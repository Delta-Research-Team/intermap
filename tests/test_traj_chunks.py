# Created by rglez at 1/27/25
import numpy as np

from intermap import commons as tt


def test_split_array():
    """Array is split into chunks of the specified size"""
    array = np.arange(10)
    chunk_size = 3
    chunks = list(tt.split_in_chunks(array, chunk_size))
    assert len(chunks) == 4
    assert np.array_equal(chunks[0], np.array([0, 1, 2]))
    assert np.array_equal(chunks[1], np.array([3, 4, 5]))
    assert np.array_equal(chunks[2], np.array([6, 7, 8]))
    assert np.array_equal(chunks[3], np.array([9]))


def test_split_empty():
    """Empty array returns no chunks"""
    array = np.array([])
    chunk_size = 3
    chunks = list(tt.split_in_chunks(array, chunk_size))
    assert len(chunks) == 0


def test_split_bigger():
    """Chunk size larger than array length returns the whole array as one chunk"""
    array = np.arange(5)
    chunk_size = 10
    chunks = list(tt.split_in_chunks(array, chunk_size))
    assert len(chunks) == 1
    assert np.array_equal(chunks[0], array)


def test_split_one():
    """Chunk size of one returns each element as a separate chunk"""
    array = np.arange(5)
    chunk_size = 1
    chunks = list(tt.split_in_chunks(array, chunk_size))
    assert len(chunks) == 5
    for i in range(5):
        assert np.array_equal(chunks[i], np.array([i]))
