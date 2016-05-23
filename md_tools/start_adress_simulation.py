#!/usr/bin/env python
"""
Copyright (C) 2016 Jakub Krajniak <jkrajniak@gmail.com>

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

import espressopp  # NOQA
import math  # NOQA
try:
    import MPI
except ImportError:
    from mpi4py import MPI
import time
import logging
import random
import shutil

import numpy as np
import gromacs_topology
import tools_sim as tools

from app_args import _args_md

# GROMACS units, kJ/mol K
kb = 0.0083144621

h5md_group = 'atoms'

__doc__ = 'Run GROMACS-like simulation'

# Do not to modify lines below.


def sort_trajectory(trj, ids):
    print('Sorting trajectory')
    idd = [
        x[1] for x in sorted([(p_id, col_id) for col_id, p_id in enumerate(ids)],
                             key=lambda y: (True, y[0]) if y[0] == -1 else (False, y[0]))
    ]
    return trj[idd]


def main():  #NOQA
    parser = _args_md()
    parser.add_argument('--ex', help='Explicity region size', type=float)
    parser.add_argument('--hy', help='Hybrid region size', type=float)
    parser.add_argument('--ex_position_x', help='Position of explicit region (x-direction)', type=float)
    args = parser.parse_args()

    parser.save_to_file('{}params.out'.format(args.output_prefix), args)

    if args.debug:
        for s in args.debug.split(','):
            print('Activating logger {}'.format(s))
            logging.getLogger(s.strip()).setLevel(logging.DEBUG)

    table_groups = map(str.strip, args.table_groups.split(','))
    lj_cutoff = args.lj_cutoff
    cg_cutoff = args.cg_cutoff
    max_cutoff = max([lj_cutoff, cg_cutoff])
    dt = args.dt

    time0 = time.time()
    input_conf = gromacs_topology.read(args.conf, args.top)

    box = (input_conf.Lx, input_conf.Ly, input_conf.Lz)
    print('Setup simulation...')

    # Tune simulation parameter according to arguments
    integrator_step = args.int_step
    sim_step = args.run / integrator_step

    if args.skin:
        skin = args.skin

    rng_seed = args.rng_seed
    if not args.rng_seed:
        rng_seed = random.randint(10, 1000000)

    print('Skin: {}'.format(skin))
    print('RNG Seed: {}'.format(rng_seed))

    part_prop, all_particles, adress_tuple = tools.genParticleList(
        input_conf, use_velocity=True, adress=True, use_charge=True)
    print('Reads {} particles with properties {}'.format(len(all_particles), part_prop))

    # Generate initial velocities, only for CG particles.
    particle_list = []
    index_adrat = part_prop.index('adrat')
    if 'v' not in part_prop:
        print('Generating velocities from Maxwell-Boltzmann distribution for T={}'.format(
            args.temperature))
        part_prop.append('v')
        cg_particles = [x for x in all_particles if x.adrat == 0]
        vx, vy, vz = espressopp.tools.velocities.gaussian(
            args.temperature, len(cg_particles), [x.mass for x in cg_particles],
            kb=kb)
        cg_id = 0
        last_vel = (0.0, 0.0, 0.0)
        for p in all_particles:
            t = list(p)
            if p.adrat == 0:
                last_vel = (vx[cg_id], vy[cg_id], vz[cg_id])
                cg_id += 1
            del t[index_adrat]
            t.append(espressopp.Real3D(last_vel))
            particle_list.append(t)
    else:
        for p in all_particles:
            t = list(p)
            del t[index_adrat]
            particle_list.append(t)

    del part_prop[index_adrat]

    density = sum(input_conf.masses)*1.6605402 / (box[0] * box[1] * box[2])
    print('Density: {} kg/m^3'.format(density))
    print('Box: {} nm'.format(box))

    system = espressopp.System()
    system.rng = espressopp.esutil.RNG(rng_seed)
    system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
    system.skin = skin
    if args.node_grid:
        nodeGrid = map(int, args.node_grid.split(','))
    else:
        nodeGrid = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)
    print('Number of nodes {}, node-grid: {}'.format(
        MPI.COMM_WORLD.size, nodeGrid))
    cellGrid = espressopp.tools.decomp.cellGrid(box, nodeGrid, max_cutoff, skin)

    print('Cell grid: {}'.format(cellGrid))

    system.storage = espressopp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

    system.storage.addParticles(map(tuple, particle_list), *part_prop)

    vs_list = espressopp.FixedVSList(system.storage)
    vs_list.addTuples(adress_tuple)

    at_particle_group = espressopp.ParticleGroup(system.storage)
    at_particle_ids = set()
    cg_particle_ids = set()
    for a in adress_tuple:
        cg_particle_ids.add(a[0])
        for at in a[1:]:
            at_particle_group.add(at)
            at_particle_ids.add(at)

    integrator = espressopp.integrator.VelocityVerletHybrid(system, vs_list)
    integrator.dt = args.dt

    system.storage.decompose()

    print('Setting AdResS simulation')
    print('Explicit region: {}, hybrid: {}'.format(args.ex, args.hy))
    print('Explicit region position: {}'.format(args.ex_position_x*box[0]))

    adress = espressopp.integrator.AdressNew(system, vs_list)
    adress.setHyEx(espressopp.Real3D(args.ex_position_x*box[0], 0, 0), args.ex, args.hy, False)
    integrator.addExtension(adress)

    exclusionlistAT = [p for p in input_conf.exclusions
                       if p[0] in at_particle_ids and p[1] in at_particle_ids]
    exclusionlistCG = [p for p in input_conf.exclusions
                       if p[0] in cg_particle_ids and p[1] in cg_particle_ids]
    print('Excluded pairs for LJ interaction (AT): {}'.format(len(exclusionlistAT)))
    print('Excluded pairs for LJ interaction (CG): {}'.format(len(exclusionlistCG)))
    verletlistAT = espressopp.VerletListHybridAT(
        system, cutoff=args.lj_cutoff, exclusionlist=exclusionlistAT)

    verletlistCG = espressopp.VerletListHybridCG(
        system, cutoff=args.cg_cutoff, exclusionlist=exclusionlistCG)

    lj_interaction = espressopp.interaction.VerletListHybridLennardJones(
        verletlistAT, False)
    lj_interaction = tools.setLennardJonesInteractions(
        system, input_conf, verletlistAT, args.lj_cutoff,
        input_conf.nonbond_params,
        interaction=lj_interaction)
    coulomb_interaction = espressopp.interaction.VerletListHybridReactionFieldGeneralized(
        verletlistAT, False)
    coulomb_interaction = gromacs_topology.setCoulombInteractions(
        system, verletlistAT, 0.9, input_conf.atomtypeparams,
        epsilon1=args.coulomb_epsilon1,
        epsilon2=args.coulomb_epsilon2, kappa=args.coulomb_kappa,
        interaction=coulomb_interaction)
    pair14_interactions = tools.setPairInteractions(
        system, input_conf, args.lj_cutoff)
    tab_cg = tools.setTabulatedInteractions(
        system, input_conf.atomtypeparams,
        vl=verletlistCG,
        cutoff=args.cg_cutoff,
        interaction=espressopp.interaction.VerletListHybridTabulated(
            verletlistCG, True))
    if lj_interaction is not None:
        system.addInteraction(lj_interaction, 'lj')
    if coulomb_interaction is not None:
        system.addInteraction(coulomb_interaction, 'coulomb')
    if tab_cg is not None:
        system.addInteraction(tab_cg, 'lj-tab')
    tools.setBondedInteractions(
        system, input_conf)
    tools.setAngleInteractions(
        system, input_conf)
    tools.setDihedralInteractions(
        system, input_conf)

    print('Bonds: {}'.format(sum(len(x) for x in input_conf.bondtypes.values())))
    print('Angles: {}'.format(sum(len(x) for x in input_conf.angletypes.values())))
    print('Dihedrals: {}'.format(sum(len(x) for x in input_conf.dihedraltypes.values())))
    print('Pairs: {}'.format(sum(len(x) for x in input_conf.pairtypes.values())))


    # Define the thermostat
    if args.temperature:
        temperature = args.temperature * kb
    else:
        temperature = 0.0
    print('Temperature: {}, gamma: {}'.format(args.temperature, args.thermostat_gamma))
    print('Thermostat: {}'.format(args.thermostat))
    thermostat = espressopp.integrator.LangevinThermostatOnGroup(system, at_particle_group)
    thermostat.temperature = temperature
    thermostat.gamma = args.thermostat_gamma
    integrator.addExtension(thermostat)

    print("Decomposing now ...")
    system.storage.decompose()

    # Observe tuple lists
    energy_file = '{}_energy_{}.csv'.format(args.output_prefix, rng_seed)
    print('Energy saved to: {}'.format(energy_file))
    system_analysis = espressopp.analysis.SystemMonitor(
        system,
        integrator,
        espressopp.analysis.SystemMonitorOutputCSV(energy_file))
    temp_comp = espressopp.analysis.Temperature(system)
    system_analysis.add_observable('T', temp_comp)
    system_analysis.add_observable(
        'Ekin', espressopp.analysis.KineticEnergy(
            system, temp_comp))
    for label, interaction in sorted(system.getAllInteractions().items()):
        print('System analysis: adding {}'.format(label))
        system_analysis.add_observable(
            label, espressopp.analysis.PotentialEnergy(system, interaction))
    ext_analysis = espressopp.integrator.ExtAnalyze(system_analysis, args.energy_collect)
    integrator.addExtension(ext_analysis)
    print('Configured system analysis')

    h5md_output_file = '{}_{}_{}'.format(args.output_prefix, rng_seed, args.output_file)
    print('Save trajectory to: {}'.format(h5md_output_file))
    traj_file = espressopp.io.DumpH5MD(
        system, h5md_output_file,
        group_name=h5md_group,
        static_box=False,
        author='Jakub Krajniak',
        email='jkrajniak@gmail.com',
        store_species=args.store_species,
        store_state=args.store_state,
        store_lambda=args.store_lambda)

    traj_file.set_parameters({
        'temperature': args.temperature})

    print('Reset total velocity')
    total_velocity = espressopp.analysis.TotalVelocity(system)
    total_velocity.reset()

    integrator.step = args.initial_step

    k_trj_collect = int(math.ceil(float(args.trj_collect) / integrator_step))
    k_energy_collect = int(math.ceil(float(args.energy_collect) / integrator_step))

    print('Running simulation for {} steps'.format(sim_step*integrator_step))
    print('Collect trajectory every {} step'.format(k_trj_collect*integrator_step))
    print('Collect energy every {} step'.format(k_energy_collect*integrator_step))

    if args.interactive:
        import IPython
        IPython.embed()

    for k in range(sim_step):
        if k_energy_collect > 0 and k % k_energy_collect == 0:
            system_analysis.info()
        if k_trj_collect > 0 and k % k_trj_collect == 0:
            int_step = args.initial_step + k*integrator_step
            traj_file.dump(int_step, int_step*dt)
        if k_trj_collect > 0 and k % 100 == 0:
            traj_file.flush()
        integrator.run(integrator_step)
    else:
        traj_file.dump(sim_step*integrator_step, sim_step*integrator_step*dt)
        traj_file.close()

    # Saves output file.
    output_gro_file = '{}_{}_confout.gro'.format(args.output_prefix, rng_seed)
    dump_gro = espressopp.io.DumpGRO(
        system, integrator, filename=output_gro_file,
        unfolded=True, append=False)
    dump_gro.dump()
    print('Wrote end configuration to: {}'.format(output_gro_file))

    print('finished!')
    print('total time: {}'.format(time.time()-time0))
    espressopp.tools.analyse.final_info(system, integrator, verletlist, time0, time.time())


if __name__ == '__main__':
    main()
