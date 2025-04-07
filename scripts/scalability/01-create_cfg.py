# Created by rglez at 4/6/25
import configparser
import os
from os.path import dirname, join

# =============================================================================
#
# =============================================================================
template = '/home/rglez/RoyHub/intermap/data/scalability_by_RGA/temp.cfg'
out_dir = dirname(template)
n_procs = [1, 2, 4, 8, 12]
chunk_sizes = [1, 100, 200, 400, 800]
resolution = ['atom', 'residue']
# =============================================================================


# Read the configuration file
config_obj = configparser.ConfigParser(allow_no_value=True,
                                       inline_comment_prefixes='#')
config_obj.optionxform = str
config_obj.read(template)

# Create the configuration files
job_names = [f'{x}-{y}-{z}'
             for x in resolution for y in n_procs for z in chunk_sizes]

with open(join(out_dir, 'runner.sh'), 'w') as f:
    f.write('#!/bin/bash\n\n')

    for job in job_names:
        # Create a copy of the template configuration
        resolution, n_procs, chunk_size = job.split('-')
        config_obj['generals']['job_name'] = job
        config_obj['generals']['n_procs'] = n_procs
        config_obj['generals']['output_dir'] = join(out_dir, job)
        config_obj['topo-traj']['chunk_size'] = chunk_size
        config_obj['interactions']['resolution'] = resolution

        # Create a directory for the job
        job_dir = config_obj['generals']['output_dir']
        os.makedirs(job_dir, exist_ok=True)

        # Write the configuration file
        with open(cfg := join(job_dir, f'{job}.cfg'), 'w') as configfile:
            config_obj.write(configfile)

        # Intermap command
        txt = cfg.replace('.cfg', '.txt')
        cmd = f"/usr/bin/time -v intermap {cfg}  2> {txt}\n"
        f.write(cmd)
