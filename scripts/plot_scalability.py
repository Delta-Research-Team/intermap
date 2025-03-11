# Created by rglez at 3/9/25
from rgpack import generals as gnl


def get_ram_time(file):
    with open(file) as f:
        lines = f.readlines()
    ram, time = 0, 0
    for line in lines:
        if line.startswith('\tUser time'):
            try:
                time = float(line.split()[3])
            except ValueError:
                time = 0
        if line.startswith('\tMaximum resident set size'):
            try:
                ram = int(line.split()[5]) / 1024
            except ValueError:
                ram = 0
    return ram, time


# =============================================================================
#
# =============================================================================
in_dir = '/media/rglez/Expansion/RoyData/intermap/benchmark/core_chunk'

# =============================================================================

chunk_files = list(gnl.recursive_finder('chunk*', in_dir))
tram = gnl.recursive_defaultdict()
for file in chunk_files:
    case = file.split('/')[-3]
    traj = file.split('/')[-2]
    sele = file.split('/')[-1]
    tram[case][traj][sele] = get_ram_time(file)
