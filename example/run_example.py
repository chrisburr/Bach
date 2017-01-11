#!/usr/bin/env python
from __future__ import division
from __future__ import print_function

import os
from os.path import join
import random
import shutil
import subprocess


MISALIGNMENTS = {
    '{Tx}': lambda: str(random.gauss(0, 0.05)),
    '{Ty}': lambda: str(random.gauss(0, 0.05)),
    '{Tz}': lambda: '0',
    '{Rx}': lambda: '0',
    '{Ry}': lambda: '0',
    '{Rz}': lambda: '0',
}


for seed in range(1, 1000):
    print('Running for seed', seed)

    out_dir = join('output', str(seed))
    os.makedirs(out_dir)
    random.seed(seed)

    # Generate a misalignment scenario
    output_lines = []
    with open('xml/conditions-template.xml') as f:
        for line in f.readlines():
            for key, get_value in MISALIGNMENTS.items():
                if key in line:
                    assert line.count(key) == 1
                    line = line.replace(key, get_value())
            output_lines.append(line)

    with open(join(out_dir, 'conditions-misaligned.xml'), 'wt') as f:
        f.writelines(output_lines)

    # Set the seed for generating the toy data
    output_lines = []
    with open('GenerateToy.xml') as f:
        for line in f.readlines():
            line = line.replace(
                '<constant name="Seed" type="I" value="42"/>',
                '<constant name="Seed" type="I" value="{}"/>'.format(seed)
            )
            output_lines.append(line)

    with open(join(out_dir, 'GenerateToy.xml'), 'wt') as f:
        f.writelines(output_lines)

    # Copy the other required files
    shutil.copy('xml/conditions-nominal.xml', out_dir)
    shutil.copy('xml/repository-nominal.xml', out_dir)
    shutil.copy('xml/repository-misaligned.xml', out_dir)
    shutil.copy('xml/default_module.xml', out_dir)
    shutil.copy('xml/Telescope.xml', out_dir)
    shutil.copy('AlignToyData.xml', out_dir)

    os.chdir(out_dir)

    # Generate the toy data
    subprocess.check_output(
        'geoPluginRun -volmgr -destroy -plugin Bach_main '
        '-config GenerateToy.xml -input file:$PWD/Telescope.xml '
        '-deltas file:$PWD/repository-nominal.xml 2>&1 | tee generate.log',
        stderr=subprocess.STDOUT, shell=True)

    # Try to align the toy data
    subprocess.check_output(
        'geoPluginRun -volmgr -destroy -plugin Bach_main '
        '-config AlignToyData.xml -input file:$PWD/Telescope.xml '
        '-deltas file:$PWD/repository-misaligned.xml 2>&1 | tee alignment.log',
        stderr=subprocess.STDOUT, shell=True)

    os.chdir('../..')
