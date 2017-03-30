#!/usr/bin/env python
"""
Copyright (C) 2017 Jakub Krajniak <jkrajniak@gmail.com>

This file is distributed under free software licence:
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
from md_libs import files_io


def _args():
    parser = argparse.ArgumentParser('Sort dihedrals into improper and proper')
    parser.add_argument('topology')

    return parser.parse_args()


def main():
    args = _args()

    in_top = files_io.GROMACSTopologyFile(args.topology)
    in_top.read()

    # improper dihedrals, func == 1
    in_top.improper_dihedrals.update({k: v for k, v in in_top.dihedrals.items() if int(v[0]) == 1})
    in_top.dihedrals = {k: v for k, v in in_top.dihedrals.items() if int(v[0]) != 1}
    in_top.write(force=True)

if __name__ == '__main__':
    main()

