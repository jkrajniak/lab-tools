import sys
import numpy as np
cimport numpy as np

cpdef inline double calc_distance_sqr(np.ndarray pos_1, np.ndarray pos_2, np.ndarray box, np.ndarray half_box):
    cdef np.ndarray d = pos_1 - pos_2
    cdef int i
    for i in range(3):
        if d[i] > half_box[i]:
            d[i] -= box[i]
        elif d[i] < -half_box[i]:
            d[i] += box[i]

    return d.dot(d)


cdef inline double calc_distance(np.ndarray pos_1, np.ndarray pos_2, np.ndarray box, np.ndarray half_box):
    return np.sqrt(calc_distance_sqr(pos_1, pos_2, box, half_box))


cdef inline np.ndarray calc_distance_vec(np.ndarray pos_1, np.ndarray pos_2, np.ndarray box, np.ndarray half_box):
    cdef np.ndarray d = pos_1 - pos_2
    for i in range(3):
        if d[i] > half_box[i]:
            d[i] -= box[i]
        elif d[i] < -half_box[i]:
            d[i] += box[i]
    return d


cdef inline np.ndarray fold_vec(np.ndarray d, np.ndarray box, np.ndarray half_box):
    for i in range(3):
        if d[i] > half_box[i]:
            d[i] -= box[i]
        elif d[i] < -half_box[i]:
            d[i] += box[i]
    return d


cpdef list calculate_bond(np.ndarray atom_tuples, np.ndarray box, np.ndarray half_box):
    cdef list bonds = []
    cdef np.ndarray d
    cdef float val
    for p in atom_tuples:
        if p[0] is not None and p[1] is not None:
            bonds.append(calc_distance(p[1], p[0], box, half_box))
        else:
            bonds.append(None)
    return bonds


cpdef list calculate_bond_vec(np.ndarray atom_tuples, np.ndarray box, np.ndarray half_box):
    cdef list bonds = []
    cdef np.ndarray d
    cdef float val
    for p in atom_tuples:
        if p[0] is not None and p[1] is not None:
            bonds.append(calc_distance_vec(p[1], p[0], box, half_box))
        else:
            bonds.append(None)
    return bonds


cpdef list calculate_bond_vec_0(np.ndarray atom_tuples, np.ndarray box, np.ndarray half_box):
    cdef list bonds = []
    cdef np.ndarray d
    cdef float val
    for p in atom_tuples:
        if p[0] is not None and p[1] is not None:
            bonds.append(p[0] - p[1])
            #bonds.append(calc_distance_vec(p[1], p[0], box, half_box))
        else:
            bonds.append(None)
    return bonds


cpdef list calculate_angle(np.ndarray atom_triplets, np.ndarray box, np.ndarray half_box):
    cdef list angles = []
    cdef double cosang
    cdef double sinang
    cdef np.ndarray p1, p2, p3, v1, v2
    for p in atom_triplets:
        if p[0] is not None and p[1] is not None and p[2] is not None:
            p1, p2, p3 = p[0], p[1], p[2]
            ## OLD
            #v1 = p2 - p1
            #v2 = p2 - p3
            #cosang = np.dot(v1, v2)
            #sinang = np.linalg.norm(np.cross(v1, v2))
            #angles.append(np.rad2deg(np.arctan2(sinang, cosang)))

            v1 = calc_distance_vec(p1, p2, box, half_box)
            v2 = calc_distance_vec(p3, p2, box, half_box)
            v1_sqr = np.dot(v1, v1)
            v2_sqr = np.dot(v2, v2)
            cos_theta = np.dot(v1, v2) / (np.sqrt(v1_sqr) * np.sqrt(v2_sqr))
            angles.append(np.rad2deg(np.arccos(cos_theta)))
            #cosang = np.dot(v1, v2)
            #sinang = np.linalg.norm(np.cross(v1, v2))
            #angles.append(np.rad2deg(np.arctan2(sinang, cosang)))
        else:
            angles.append(None)
    return angles


cpdef list calculate_dihedral(np.ndarray atom_quadruples, np.ndarray box, np.ndarray half_box):
    cdef list dihedrals = []
    cdef np.ndarray b, v, b1, m, r21, r32, r43, rijk, rjkn
    cdef double x, y, rijk_sqr, rjkn_sqr, rijk_abs, rjkn_abs, inv_rijk, inv_rjkn, cos_phi
    for p in atom_quadruples:
        if p[0] is not None and p[1] is not None and p[2] is not None and p[3] is not None:
            b = np.array([calc_distance_vec(p[i], p[i+1], box, half_box) for i in range(3)])
            b[0] *= -1
            v = np.array([v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]]])
            v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1, 1)
            b1 = b[1] / np.linalg.norm(b[1])
            x = np.dot(v[0], v[1])
            m = np.cross(v[0], b1)
            y = np.dot(m, v[1])
            dihedrals.append(np.rad2deg(np.arctan2(y, x)))
        else:
            dihedrals.append(None)
    return dihedrals


cpdef process_distance(np.ndarray atom_pairs, np.ndarray box, np.ndarray half_box):
    bonds = []
    cdef int l = len(atom_pairs)
    cdef Py_ssize_t i
    for i in range(l):
        pair = atom_pairs[i]
        if len(pair) == 2:
            bonds.append(calc_distance(pair[0], pair[1], box, half_box))
    return bonds


cpdef tuple calculate_end_end_acf(np.ndarray trajectory, int max_tau):
    """Calculating autocorrelation function of end-to-end vector
    Args:
        trajectory: numpy array with trajectory.
        max_tau: Maximum tau.
    Returns:
        Returns ACF and ACF error.
    """
    cdef int N = len(trajectory[0])  # number of molecules
    cdef int max_t = len(trajectory)
    cdef int m, n, i

    cdef np.ndarray acf = np.zeros(max_tau)
    cdef np.ndarray acf_error = np.zeros(max_tau)
    cdef double sumdist
    cdef list intermediate_results
    print 'Frames: 0..', max_t
    print 'Max tau:', max_tau
    # Go over time.
    for m in xrange(max_tau):
        sys.stdout.write('Tau: {}\r'.format(m))
        sys.stdout.flush()
        intermediate_results = []
        for n in xrange(0, max_t - m):
             sumdist = 0.0
             for i in xrange(N):
                 a1 = trajectory[n][i]
                 a2 = trajectory[n+m][i]
                 sumdist += np.dot(a1, a2)
             intermediate_results.append(sumdist / N)
        acf[m] = np.average(intermediate_results)
        acf_error[m] = np.std(intermediate_results, ddof=1) #/np.sqrt(max_t-m)
    return acf, acf_error


cpdef tuple calculate_end_end_acf_single(np.ndarray trajectory, int tau):
    """Calculating autocorrelation function of end-to-end vector
    Args:
        trajectory: numpy array with trajectory.
        tau: tau.
    Returns:
        Returns ACF and ACF error.
    """
    cdef int N = len(trajectory[0])  # number of molecules
    cdef int max_t = len(trajectory)
    cdef int m, n, i

    cdef double acf = 0.0
    cdef double acf_error = 0.0
    cdef double sumdist = 0.0
    cdef list intermediate_results = []
    
    # Go over time.
    sys.stdout.write('Tau: {}\r'.format(tau))
    sys.stdout.flush()
    for n in xrange(0, max_t - m):
         sumdist = 0.0
         for i in xrange(N):
             sumdist += np.dot(trajectory[n][i], trajectory[n+m][i])
         intermediate_results.append(sumdist / N)
    acf = np.average(intermediate_results)
    acf_error = np.std(intermediate_results, ddof=1) #/np.sqrt(max_t-m)
    return acf, acf_error


cpdef list calculate_msd_internal_distance(np.ndarray frame, int nchains, int chain_length, np.ndarray box, np.ndarray half_box):
    """Calculates mean square internal distance.
    Args:
        frame: single time frame.
        nchains: Number of polymer chains
        chain_length: Length of polymer chain
        box: The box.
        half_box: The half size of box.
    Returns:
        The data set indexed by n.
    """
    cdef np.ndarray output = np.zeros(chain_length)
    cdef np.ndarray error = np.zeros(chain_length)
    cdef int n, j, t, i
    cdef int m
    cdef double sumsqdist = 0.0
    cdef double intdist_sum = 0.0
    cdef np.ndarray d = np.zeros(3)
    cdef list intermediate_results = [[] for _ in range(chain_length)]

    for n in xrange(0, chain_length):
        intdist_sum = 0.0
        for j in xrange(nchains):
            sumsqdist = 0.0
            for i in xrange(j*chain_length, j*chain_length + chain_length - n):
                d = frame[i+n] - frame[i]
                sumsqdist += d.dot(d)
            intdist_sum += sumsqdist  / (chain_length - n)
        intermediate_results[n].append(intdist_sum/nchains)
    return intermediate_results


cpdef tuple calculate_com_chains(np.ndarray traj, int chain_length, int chains, np.ndarray masses, double tot_mass):
    """Calculates COM of chains.
    Args:
        traj: single time frame.
        chain_length: length of chain.
        chains: Number of chains.
        masses: The numpy array with mass for particles in single chain.
        tot_mass: The total mass of single chain
    Returns:
        The tuple with center of mass of chains and center of mass of whole system.
    """
    cdef np.ndarray output = np.zeros(shape=(chains, 3))
    cdef np.ndarray output_sys = np.zeros(3)
    cdef np.ndarray com = np.zeros(3)
    cdef np.ndarray sys_com = np.zeros(3)

    for ch in xrange(chains):
        output[ch] = np.zeros(3)
        for n in xrange(chain_length):
            output[ch] += traj[ch*chain_length+n]*masses[n]
        output[ch] /= tot_mass
        sys_com += output[ch]
    output_sys = sys_com / chains
    return output, output_sys


cpdef tuple calculate_msd_single(np.ndarray trj_com, np.ndarray trj_sys_com, np.ndarray box, int N, int r_every, int tau):
    """Calculate MSD for given value of tau (m).
    Args:
        trj_com: numpy array with trajectory of com.
        trj_sys_com: numpy array with trajectory of system com.
        box: box size.
        N: number of chains.
        tau: Given tau.
        r_every: restart every .. time frame.
    Returns:
        MSD and MSD error
    """
    cdef int n, j, i
    cdef double sumdist
    cdef int max_t = len(trj_com)
    
    #cdef np.ndarray inv_box = 1.0/box

    sys.stdout.write('tau: %d\n' % tau)

    cdef np.ndarray pos1, pos2
    cdef list intermediate_results

    ror = ['-', '\\', '|', '/', '-', '\\', '|', '/', '-', '\\']

    intermediate_results = []
    for n in xrange(0, max_t - tau, r_every):  # starting from t=n ... max_t - tau
        sumdist = 0.0
        for i in xrange(N):
            pos1 = trj_com[n+tau][i] - trj_sys_com[n+tau]
            pos2 = trj_com[n][i] - trj_sys_com[n]
            d = pos2 - pos1
            #d -= np.round(d*inv_box)*box
            sumdist += d.dot(d)
        intermediate_results.append(sumdist / N)
        #sys.stdout.write('%s\r' % ror[n % 10])
        #sys.stdout.flush()
    #sys.stdout.write('\n')
    return np.average(intermediate_results), np.std(intermediate_results, ddof=1)
