# SDPD-mf_thixo

<img width="1023" height="670" alt="image" src="https://github.com/user-attachments/assets/2231c7d1-8886-4bb9-adc5-8892ea620e60" />

Smoothed Dissipative Particle Dynamics (SDPD) model for multiphase and thixotropic fluids implemented in LAMMPS.

# Description

An extended version of the standard SDPD model in LAMMPS (Jalalvand et al. 2021),  including:
 	
a) Update non-slip in flat wall: This update included the option for a non-slip boundary condition in a flat wall as described in Bian and Ellero 2012. The correction can be turned on and off using the variable slip[k][l], where the index corresponds to the type of particles k and l. Thus, a fluid particle type k can interact with a non-slip bc with a wall particle type l.

b) Update MF: This updated multiphase (mf) version is based on the description presented in Lei et al. 2016 (10.1103/PhysRevE.94.023304), including a pair-force contribution between particles.

c) Transient viscosity model: This update includes a thixotropic (viscosity transient) model to simulate complex multiphase flows

# Installation

1. Replace the files atom.cpp, atom.h, atom_vec.cpp, atom_vec.h, fix_sph.cpp and fix_sph.h with the files from the .tar
2. Copy the files pair_sdpd_taitwater_isothermal_mf.cpp and pair_sdpd_taitwater_isothermal_mf.h in the folder DPD-SMOOTH
3. Run: make mpi

# Examples

Here you will find an example (folder named “Example”) containing the files needed to run a case of the interaction of three droplets under a continuous fluid phase. The main script is called “in.sdpd_phase.2d”.
