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
import math
import numpy as np
import sys

from matplotlib import pyplot as plt

col = int(sys.argv[2])
input_data = np.loadtxt(sys.argv[1])[:, col]*(-0.101325)

n = input_data.shape[0]
out = []

for tb in range(4, 30000, 100):
    nb = int(math.ceil(n/tb))
    #print tb, nb, input_data.shape[0], tb*nb
    data = input_data[:tb*nb]
    nblocks = np.split(data, nb)
    mean = np.mean(data)
    tot_var = np.var(data)
    avg_block = np.average(nblocks, axis=1)
    var_block = np.mean(np.power(avg_block - mean, 2))
    out.append([tb, tb*var_block/tot_var])

out = np.array(out)
plt.plot(out[:, 0], out[:, 1])
plt.show()
