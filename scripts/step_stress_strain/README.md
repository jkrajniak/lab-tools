Stress-strain calculation with static method
==============================================

This is a semi-automate way of calculating stress-strain curve
with the static method.

The `stress_strain_step.sh` file is prepared to extend the
box in *z-direction* by 1% every step, then perform NPT simulation
to collect data at that strain point.

The `make_deform.sh` takes the initial box from `conf.gro` file
and return .mdp file where the deformation of the box in *z-direction*
is made during 500ps of NPT simulation.
