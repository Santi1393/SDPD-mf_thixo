# SDPD-mf_thixo

<img width="1023" height="670" alt="image" src="https://github.com/user-attachments/assets/2231c7d1-8886-4bb9-adc5-8892ea620e60" />

Smoothed Dissipative Particle Dynamics (SDPD) model for multiphase and thixotropic fluids implemented in LAMMPS.

# Description

An extended version of the standard SDPD model in LAMMPS (Jalalvand et al. 2021),  including:
 	
a) Update non-slip in flat wall: This update included the option for a non-slip boundary condition in a flat wall as described in Bian and Ellero 2012. The correction can be turned on and off using the variable slip[k][l], where the index corresponds to the type of particles k and l. Thus, a fluid particle type k can interact with a non-slip bc with a wall particle type l.

b) Update MF: This updated multiphase (mf) version is based on the description presented in Lei et al. 2016 (10.1103/PhysRevE.94.023304), including a pair-force contribution between particles.

c) Transient viscosity model: This update includes a thixotropic (viscosity transient) model to simulate complex multiphase flows

# Software version

This SDPD extension must be compiled in the 21-Nov-2023 release of the open source LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator).

# Installation and Build flags

The installation process and build flags are detailed step by step below.

1. Download and unzip the .tar.gz file “lammps-21Nov23.tar.gz” found in this repository.
2. Replace the following files with the files from the .tar

   *atom.cpp*
   
   *atom.h*
   
   *atom_vec.cpp*
   
   *atom_vec.h*
   
   *fix_sph.cpp*
   
   *fix_sph.h*

4. Copy the files *pair_sdpd_taitwater_isothermal_mf.cpp* and *pair_sdpd_taitwater_isothermal_mf.h* in the folder "DPD-SMOOTH"
5. In the same folder, locate the “Install.sh” file. Include these two lines to indicate that you are going to install a new dependency:

   *action pair_sdpd_taitwater_isothermal_mf.h*
   
   *action pair_sdpd_taitwater_isothermal_mf.cpp*
   
6. LAMMPS must be compiled with the SPH and SDPD package enabled. To activate it, located in the “src” folder type:

   *make yes-SPH*
   
   *make yes-DPD-SMOOTH*
 
8. LAMMPS must be compiled using MPI support. Located in the “src” folder type:

   *make pu*
   
   *make mpi*

10. Run a numerical example to verify the installation. Copy any example to your local folder and, once located in the case folder, type:

   *lmp_mpi <in.sdpd_phase.2d*

# Numerical examples

In the folder called “Numerical examples” you will find a collection of all the numerical cases covered in this research. Detailed information on each case and how to run them can be found in the README file in each example folder. The complete list of cases studied is detailed below: 

1. Surface tension of a droplet - Used for static validation - See Figure 1_(a)
2. Retraction of a stretched droplet - Used for static validation - See Figure 1_(b)
3. Static contact angle between droplet and solid wall - Used for static validation - See Figure 1_(c)
4. Triple contact angle between three diferents droplets - Used for static validation - See detail in Figure 1_(c)
5. Poiseuille flow for one phase - Used for dynamic validation - See Figure A.12_(a) and A.12_(d)
6. Reverse Poiseuille flow for one phase - Used for dynamic validation - See Figure A.12_(b) and A.12_(e)
7. Flow around a cylinder - Used for dynamic validation - See Figure A.12_(c) and A.12_(f)
8. Poiseuille flow for two phases - Used for dynamic validation - See Figure 2_(a)
9. Taylor deformation vs capillary number - Used for dynamic validation - See Figure 2_(b)
10. Droplet break-up - Used for dynamic validation - See Figure 2_(c)
11. Thixotropic shear flow in a channel - Used for thixotropic validation - See Figure 3_(a) and 3_(b)
12. Liquid-Liquid Phase Separation (LLPS) - First exploratory case - See Figure 4 and 5
13. Poiseuille flow for two phases with one thixotropic phase - Second exploratory case - See Figure 6
14. Emulsion with thixotropic continuous phase and newtonian droplets - Second exploratory case - See Figure 7
15. Emulsion with newtonian continuous phase and thixotropic droplets - Second exploratory case - See Figure 8
16. Droplet dynamics in a periodically constricted channel - Third exploratory case - See Figure 9 and 10
17. Droplet merging using micro-devices - fourth exploratory case - See Figure 11



