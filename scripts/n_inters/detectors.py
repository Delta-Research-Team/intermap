# Created by rglez at 4/18/25
import os

from bitarray import bitarray as ba
from rgpack.generals import recursive_defaultdict as rd


class Detector:
    """
    Detector class to detect n_inters and n_pairs from a csv file.
    """

    def __init__(self, inters_file, ignores):
        self.n_pairs = 0
        self.n_inters = 0
        self.inters_file = inters_file
        self.ignores = ignores
        self.dejavu = set()

    def detect(self):
        """
        Detect intermap n_inters and n_pairs from a csv file.
        """
        with open(self.inters_file) as f:
            for line in f:
                start_bad = line.startswith(self.ignores)
                if start_bad:
                    continue
                else:
                    self.detect_inters(line)
        return self.n_inters, self.n_pairs

    def detect_inters(self, line):
        """
        Detect n_inters and n_pairs from a file.

        Args:
            line: string with the line to be processed
        """
        raise NotImplementedError(
            'This method should be implemented in a subclass')


class IMapDetector(Detector):
    """
    IMapDetector class to detect n_inters and n_pairs from a csv file.
    """

    def __init__(self, inters_file, ignores):
        super().__init__(inters_file, ignores)

    def detect_inters(self, line):
        """
        Detect n_inters and n_pairs from a file.

        Args:
            line: string with the line to be processed
        """
        time = line.split(',')[-1]
        bitarr = ba(time)
        self.n_inters += bitarr.count()
        self.n_pairs += 1


class GetContactsDetector(Detector):
    """
    GetContactsDetector class to detect n_inters and n_pairs from a csv file.
    """

    def __init__(self, inters_file, ignores):
        super().__init__(inters_file, ignores)

    def detect_inters(self, line):
        """
        Detect n_inters and n_pairs from a file.

        Args:
            line: string with the line to be processed
        """
        pair_line = ''.join([x for x in line.split("\t")[2:4]]).strip()
        self.dejavu.add(pair_line)
        self.n_inters += 1
        self.n_pairs = len(self.dejavu)


class ProlifDetector(Detector):
    """
    ProlifDetector class to detect n_inters and n_pairs from a file.
    """

    def __init__(self, inters_file, ignores):
        super().__init__(inters_file, ignores)

    def detect_inters(self, line):
        """
        Detect n_inters and n_pairs from a file.

        Args:
            line: string with the line to be processed
        """
        time = [0 if x == 'False' else 1 for x in line.split(',')[3:]]
        bitarr = ba(time)
        self.n_inters += bitarr.count()
        self.n_pairs += 1


class MDLRDetector(Detector):
    """
    MDLRDetector class to detect n_inters and n_pairs from a file.
    """

    def __init__(self, inters_file, ignores):
        super().__init__(inters_file, ignores)

    def detect_inters(self, line):
        """
        Detect n_inters and n_pairs from a file.

        Args:
            line: string with the line to be processed
        """
        line_info = line.split(',')
        pair = f'{line_info[0]}{line_info[1]}'
        self.dejavu.add(pair)
        time = [0 if x == '0' else 1 for x in line_info[2:]]
        bitarr = ba(time)
        self.n_inters += bitarr.count()
        self.n_pairs = len(self.dejavu)


# =============================================================================
# files
# =============================================================================
user = os.getenv('USER')
detectors = {'plif': ProlifDetector,
             'imap': IMapDetector,
             'gc': GetContactsDetector,
             'mdlr': MDLRDetector}

ignores = {'plif': ('ligand',),
           'imap': ('#', 'sel1'),
           'gc': ('#',),
           'mdlr': ('A', 'U', 'T', 'I', 'L', 'T', ' ', ',', '\n')}

traj_results = {
    'nucleosome': {
        'plif': f'/media/gonzalezroy/Roy2TB/RoyData/intermap/benchmark/n_inters/prolif/ncp-OCS-nosalt/nucleic-protein/interacciones.csv',
        # 'imap': f'/home/{user}/RoyHub/intermap/tests/imaps/prot-nuc/8oxoGA2_1_InterMap_full.csv',
        'gc': f'/media/{user}/Roy2TB/RoyData/intermap/benchmark/n_inters/getContacts/ncp-OCS-nosalt/nucleic-protein/results.tsv',
        'mdlr': None},

    'T2': {
        'plif': None,
        'imap': None,
        'gc': None,
        'mdlr': None},

    'T3': {
        'plif': None,
        'imap': None,
        'gc': None,
        'mdlr': None},

    'T4': {
        'plif': None,
        'imap': None,
        'gc': None,
        'mdlr': None},
}

# =============================================================================
# Flow
# =============================================================================
counts = rd()
for traj_label in traj_results.keys():
    for soft_case in traj_results[traj_label].keys():
        if traj_results[traj_label][soft_case] is not None:
            detector = detectors[soft_case]
            ignores_ = ignores[soft_case]
            inters_file = traj_results[traj_label][soft_case]

            self = detector(inters_file, ignores_)
            n_inters, n_pairs = self.detect()
            counts[traj_label][soft_case]['n_inters'] = n_inters
            counts[traj_label][soft_case]['n_pairs'] = n_pairs

# =============================================================================
# imap
# =============================================================================
# imap_csv = '/home/rglez/RoyHub/intermap/scripts/scalability/atom-2-24/atom-2-24_InterMap_full.csv'
# self = IMapDetector(imap_csv, ('#', 'sel1'))
# n_inters, n_pairs = self.detect()

# =============================================================================
# getContacts
# =============================================================================
gc_tsv = '/media/rglez/Roy2TB/RoyData/intermap/N_INTERS-FINAL/p53/results.tsv'
self = GetContactsDetector(gc_tsv, ('#',))
n_inters, n_pairs = self.detect()

# =============================================================================
# Prolif
# =============================================================================
# plif_csv = '/home/rglez/RoyHub/intermap/data/benchmark/n_inters/prolif/1kx5sno_dry/protein-protein/interacciones.csv'
# self = ProlifDetector(plif_csv, ('ligand',))
# n_inters, n_pairs = self.detect()

# =============================================================================
# MDLR
# =============================================================================
# mdlr_csv = '/home/rglez/RoyHub/intermap/scripts/n_inters/IFP_results_rglez_20250306_101320.csv'
# self = MDLRDetector(mdlr_csv, ('A', 'U', 'T', 'I', 'L', 'T', ' ', ',', '\n'))
# n_inters, n_pairs = self.detect()

# from bitarray import util as bu
# from rgpack import generals as gnl
# pkl = '/media/gonzalezroy/Roy2TB/RoyData/intermap/benchmark/n_inters/intermap/ncp-OCS-nosalt/nucleic-protein/nucleic-protein_InterMap.pickle'
# dicty = gnl.unpickle_from_file(pkl)
# counter = 0
# for x in dicty:
#     bitarr = bu.sc_decode(dicty[x]['time'])
#     counter += bitarr.count()


a = ba(10)
a[:5][0, 2, 4] = 1
