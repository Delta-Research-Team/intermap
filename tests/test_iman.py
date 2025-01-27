# Created by gonzalezroy at 1/24/25

import pytest
from rgpack import generals as gnl

from intermap import indices as idx


@pytest.fixture(scope="module")
def imans(topo_trajs_fix):
    """
    Prepare IndexManager objects for the tests
    """
    # Atomic selection definitions
    seles_all = {'wat': 'resname WAT',
                 'lig': 'resname LIG',
                 'prot': 'protein',
                 'nuc': 'nucleic',
                 'helix': 'resnum 204:249',
                 'ball': 'resnum 19:110',
                 'ions': 'resnum 1244:1796',
                 }

    # Pairs of selections for each trajectory
    seles_cases = {
        'trj1': [('lig', 'lig'), ('lig', 'prot'), ('lig', 'helix'),
                 ('helix', 'helix')],
        'trj2': [('ball', 'ball')],
        'trj3': [('nuc', 'wat'), ('ions', 'nuc')]}

    imans = gnl.recursive_defaultdict()
    for case, seles in seles_cases.items():
        topo, traj = topo_trajs_fix[case]
        for sel1_raw, sel2_raw in seles:
            sel1 = seles_all[sel1_raw]
            sel2 = seles_all[sel2_raw]
            label = f'{case}:{sel1_raw}:{sel2_raw}:all'
            iman = idx.IndexManager(topo, traj, sel1, sel2, 'all')
            imans[label] = iman
    return imans


def test_iman_smarts(imans):
    """
    """
    expected_values = {
        'trj1:lig:lig:all':
            {'sel_idx': 79, 'sel1_idx': 79, 'sel2_idx': 79, 'anions': 0,
             'cations': 1, 'hydroph': 22, 'metal_acc': 5, 'metal_don': 0,
             'radii': 79, 'rings': 3, 'hb_A': 5, 'hb_D': 4, 'hb_H': 4,
             'xb_A': 46, 'xb_D': 0, 'xb_H': 0},

        'trj1:lig:prot:all':
            {'sel_idx': 5067, 'sel1_idx': 79, 'sel2_idx': 4988, 'anions': 51,
             'cations': 59, 'hydroph': 914, 'metal_acc': 726, 'metal_don': 0,
             'radii': 5067, 'rings': 57, 'hb_A': 722, 'hb_D': 227, 'hb_H': 227,
             'xb_A': 2530, 'xb_D': 0, 'xb_H': 0},

        'trj1:lig:helix:all':
            {'sel_idx': 854, 'sel1_idx': 79, 'sel2_idx': 775, 'anions': 6,
             'cations': 15, 'hydroph': 156, 'metal_acc': 116, 'metal_don': 0,
             'radii': 854, 'rings': 11, 'hb_A': 116, 'hb_D': 48, 'hb_H': 48,
             'xb_A': 442, 'xb_D': 0, 'xb_H': 0},

        'trj1:helix:helix:all':
            {'sel_idx': 775, 'sel1_idx': 775, 'sel2_idx': 775, 'anions': 6,
             'cations': 14, 'hydroph': 134, 'metal_acc': 111, 'metal_don': 0,
             'radii': 775, 'rings': 8, 'hb_A': 111, 'hb_D': 44, 'hb_H': 44,
             'xb_A': 396, 'xb_D': 0, 'xb_H': 0},

        'trj2:ball:ball:all':
            {'sel_idx': 1436, 'sel1_idx': 1436, 'sel2_idx': 1436, 'anions': 28,
             'cations': 25, 'hydroph': 216, 'metal_acc': 224, 'metal_don': 0,
             'radii': 1436, 'rings': 7, 'hb_A': 224, 'hb_D': 76, 'hb_H': 76,
             'xb_A': 693, 'xb_D': 0, 'xb_H': 0},

        'trj3:nuc:wat:all':
            {'sel_idx': 448243, 'sel1_idx': 9292, 'sel2_idx': 438951,
             'anions': 873, 'cations': 0, 'hydroph': 1676, 'metal_acc': 148511,
             'metal_don': 0, 'radii': 448243, 'rings': 440, 'hb_A': 148511,
             'hb_D': 293242, 'hb_H': 293242, 'xb_A': 445707, 'xb_D': 0,
             'xb_H': 0},

        'trj3:ions:nuc:all':
            {'sel_idx': 9845, 'sel1_idx': 553, 'sel2_idx': 9292,
             'anions': 1075, 'cations': 351, 'hydroph': 1676,
             'metal_acc': 2396, 'metal_don': 0, 'radii': 9845, 'rings': 440,
             'hb_A': 2194, 'hb_D': 608, 'hb_H': 608, 'xb_A': 6756, 'xb_D': 0,
             'xb_H': 0}
    }

    for label in imans:
        assert len(imans[label].sel_idx) == expected_values[label]['sel_idx']
        assert len(imans[label].sel1_idx) == expected_values[label]['sel1_idx']
        assert len(imans[label].sel2_idx) == expected_values[label]['sel2_idx']
        assert len(imans[label].anions) == expected_values[label]['anions']
        assert len(imans[label].cations) == expected_values[label]['cations']
        assert len(imans[label].hydroph) == expected_values[label]['hydroph']
        assert len(imans[label].metal_acc) == expected_values[label]['metal_acc']
        assert len(imans[label].metal_don) == expected_values[label]['metal_don']
        assert len(imans[label].radii) == expected_values[label]['radii']
        assert len(imans[label].rings) == expected_values[label]['rings']
        assert len(imans[label].hb_A) == expected_values[label]['hb_A']
        assert len(imans[label].hb_D) == expected_values[label]['hb_D']
        assert len(imans[label].hb_H) == expected_values[label]['hb_H']
        assert len(imans[label].xb_A) == expected_values[label]['xb_A']
        assert len(imans[label].xb_D) == expected_values[label]['xb_D']
        assert len(imans[label].xb_H) == expected_values[label]['xb_H']
