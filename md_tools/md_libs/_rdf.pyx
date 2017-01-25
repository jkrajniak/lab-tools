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

Partially this file is stolen from https://github.com/pdebuyl/md_tools
and modified but the original file has the following copyright header:

 Copyright 2014 Pierre de Buyl

 This file is part of md_tools

 md_tools is free software and is licensed under the modified BSD license (see
 LICENSE file).
"""

import numpy as np
cimport numpy as np
import itertools
from libc.math cimport sqrt, floor
cimport cython

def compute_rdf(r1, r2, L, N, cutoff, do_normalize=False):
    """
    Compute the radial distribution function (rdf).

    Arguments
    ---------

    r1: [num,3] array of positions of first group
    r2: [num, 3] array of positions of second group
    L: [3] sides of a cuboid box
    N: nbins
    do_normalize: Do normalization ?

    Returns
    -------

    dx is the radius step
    all_rdf is a [N_rdf,N]
    """
    cdef int i, n_rdf, n_idx,npart
    cdef double[3] cy_L
    cdef double vol,norm,phi
    vol = 1.0
    for i in range(3):
        cy_L[i] = L[i]
        vol *= L[i]

    if cutoff == -1:
        x_max = np.min(L)/2.
    else:
        x_max = cutoff
    dx = x_max/N

    if r2 is not None:
        result = _compute_rdf_multi(r1, r2, cy_L, N, dx)
    else:
        result = _compute_rdf(r1, cy_L, N, dx)

    result = np.asarray(result)
    #if do_normalize:
    #    phi = npart/vol
    #    norm = dx*phi*npart
    #    result /= norm
    return dx, result

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef _compute_rdf_multi(double[:, ::1] r1, double[:, ::1] r2, double[3] L, int N, double dx):
    cdef int i, j, coord, idx, si, sj, rdf_idx, n_idx
    cdef double dist
    cdef double dist_sqr
    cdef double inv_dx = 1./dx
    cdef double x_max_sqr = (N*dx)**2
    cdef double pi = np.pi
    cdef double k

    cdef double[::1] result = np.zeros(N)

    for i in range(r1.shape[0]):
        si = 0
        for j in range(r2.shape[0]):
            dist_sqr = 0.0
            if all([r1[i, c] == r2[j, c] for c in range(3)]):
                continue
            for coord in range(3):
                dist = r1[i,coord]-r2[j,coord]
                if dist<-L[coord]/2.:
                    dist += L[coord]
                elif dist>L[coord]/2.:
                    dist -= L[coord]
                dist_sqr += dist**2
            if dist_sqr <= x_max_sqr:
                idx = int(floor(sqrt(dist_sqr)*inv_dx))
                result[idx] += 1

    k = 4.0*pi
    for i in range(result.shape[0]):
        result[i] /= (k*((i+0.5)*dx)**2)

    return result


def compute_rdf_index(pos, index_pairs, L, N, cutoff, do_normalize=False):
    cdef int i, n_rdf, n_idx,npart
    cdef double[3] cy_L
    cdef double vol,norm,phi
    for i in range(3):
        cy_L[i] = L[i]

    if cutoff == -1:
        x_max = np.min(L)/2.
    else:
        x_max = cutoff
    dx = x_max/N

    result = _compute_rdf_multi_idx(pos, index_pairs, cy_L, N, dx)

    result = np.asarray(result)
    #if do_normalize:
    #    phi = npart/vol
    #    norm = dx*phi*npart
    #    result /= norm
    return dx, result


cpdef compute_pairs(int[::1] input_pairs):
    return np.asarray([x for l in list(itertools.combinations_with_replacement(input_pairs, 2)) for x in l], dtype=np.int32)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef _compute_rdf_multi_idx(double[:, ::1] pos, int[::1] index_pairs, double[3] L, int N, double dx):
    cdef int i, j, coord, idx, si, sj, rdf_idx, n_idx
    cdef double dist
    cdef double dist_sqr
    cdef double inv_dx = 1./dx
    cdef double x_max_sqr = (N*dx)**2
    cdef double pi = np.pi
    cdef double k

    cdef double[3] L2
    for i in range(3):
        L2[i] = 0.5*L[i]

    cdef double[::1] result = np.zeros(N)

    cdef int pidx
    for pidx in range(index_pairs.shape[0]/2):
        i = index_pairs[pidx]
        j = index_pairs[pidx+1]
        dist_sqr = 0.0
        for coord in range(3):
            dist = pos[i, coord] - pos[j, coord]
            if dist < -L2[coord]:
                dist += L[coord]
            elif dist > L2[coord]:
                dist -= L[coord]
            dist_sqr += dist**2
        if dist_sqr <= x_max_sqr:
            idx = int(floor(sqrt(dist_sqr)*inv_dx))
            result[idx] += 1

    k = 4.0*pi
    for i in range(result.shape[0]):
        result[i] /= (k*((i+0.5)*dx)**2)

    return result


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef _compute_rdf(double[:, ::1] r1, double[3] L, int N, double dx):
    cdef int i, j, coord, idx, si, sj, rdf_idx, n_idx
    cdef double dist
    cdef double dist_sqr
    cdef double inv_dx = 1./dx
    cdef double x_max_sqr = (N*dx)**2
    cdef double pi = np.pi
    cdef double k

    cdef double[::1] result = np.zeros(N)

    for i in range(r1.shape[0]):
        for j in range(i+1, r1.shape[0]):
            dist_sqr = 0.0
            for coord in range(3):
                dist = r1[i,coord]-r1[j,coord]
                if dist<-L[coord]/2.:
                    dist += L[coord]
                elif dist>L[coord]/2.:
                    dist -= L[coord]
                dist_sqr += dist**2
            if dist_sqr <= x_max_sqr:
                idx = int(floor(sqrt(dist_sqr)*inv_dx))
                result[idx] += 1

    k = 2.0*pi
    for i in range(result.shape[0]):
        result[i] /= (k*((i+0.5)*dx)**2)

    return result

cdef inline distance_sqr(double[::1] d1, double[::1] d2, double[3] L, double[3] L2):
    cdef double dist = 0.
    cdef double dist_sqr = 0.
    cdef int coord
    for coord in range(3):
        dist = d1[coord] - d2[coord]
        if dist < -L2[coord]:
            dist += L[coord]
        elif dist > L2[coord]:
            dist -= L[coord]
        dist_sqr += dist**2
    return dist_sqr



def compute_nb(r1, r2, L, cutoff):
    cdef int i
    cdef double[3] cy_L
    cdef double[3] cy_L2
    for i in range(3):
        cy_L[i] = L[i]
        cy_L2[i] = 0.5*L[i]

    return _compute_nb(r1, r2, cy_L, cy_L2, cutoff**2)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef _compute_nb(double[:, ::1] r1, double[:, ::1] r2, double[3] L, double[3] L2, double cutoff_sqr):
    cdef int i, j, coord, idx
    cdef double dist_sqr

    cdef int[::1] result = np.zeros(r1.shape[0], dtype=np.int32)

    for i in range(r1.shape[0]):
        for j in range(r2.shape[0]):
            dist_sqr = 0.0
            for coord in range(3):
                dist = r1[i,coord]-r2[j,coord]
                if dist<-L2[coord]:
                    dist += L[coord]
                elif dist>L2[coord]:
                    dist -= L[coord]
                dist_sqr += dist**2
            if dist_sqr > 0.0 and dist_sqr <= cutoff_sqr:
                result[i] += 1

    return result

