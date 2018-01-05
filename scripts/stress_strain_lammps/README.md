Measurement of stress-strain relation with LAMMPS
==================================================

Structure 

 - `eq/` - this contains the files required to equilibrated the sample
 - `pull_step_x|y|z` - pulling in three directions (equilibrated mode)
 - `pull_x|y|z` - pulling in three directions (continuous mode)


Equilibrated mode
-------------------

The whole deformation range is divided into 50 steps. Each of the steps consists of the deformation, equilibration at given strain and measure of the stress in given direction.

Continuous mode
------------------

In this type of measurment, the sample is continuously deformed in given direction and simultaneously the stress is measured.

