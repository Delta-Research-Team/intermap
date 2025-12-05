# Created by rglez at 4/6/25
import configparser
import os
from os.path import join

from run_prolif import selections, topo_trajs

# =============================================================================
#
# =============================================================================
template = '/home/rglez/RoyHub/intermap/scripts/identity/imap.cfg'
out_dir = '/home/rglez/RoyHub/intermap/data/identity/runs/imap'
os.makedirs(out_dir, exist_ok=True)
resolutions = ['atom', 'residue']
# =============================================================================

# Read the configuration file
config_obj = configparser.ConfigParser(allow_no_value=True,
                                       inline_comment_prefixes='#')
config_obj.optionxform = str
config_obj.read(template)

with open(join(out_dir, 'runner.sh'), 'w') as f:
    f.write('#!/bin/bash\n\n')

    for case in selections:
        for resolution in resolutions:
            job_name = f"{case.split('_')[0]}-{resolution}"
            job_dir = join(out_dir, job_name)
            topology = topo_trajs[case]['topo']
            trajectory = topo_trajs[case]['traj']
            sele1 = selections[case]['sel1']
            sele2 = selections[case]['sel2']

            # Create a copy of the template configuration
            config_obj['generals']['job_name'] = job_name
            config_obj['generals']['output_dir'] = job_dir
            config_obj['interactions']['resolution'] = resolution
            config_obj['interactions']['selection_1'] = sele1
            config_obj['interactions']['selection_2'] = sele2
            config_obj['topo-traj']['topology'] = topology
            config_obj['topo-traj']['trajectory'] = trajectory

            # Create a directory for the job
            os.makedirs(job_dir, exist_ok=True)

            # Write the configuration file
            with open(cfg := join(job_dir, f'{job_name}.cfg'),
                      'w') as configfile:
                config_obj.write(configfile)

            # Intermap command
            txt = cfg.replace('.cfg', '.txt')
            cmd = f"/usr/bin/time -v intermap {cfg}  2> {txt}\n"
            f.write(cmd)
