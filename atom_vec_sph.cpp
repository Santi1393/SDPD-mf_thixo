/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 https://www.lammps.org/, Sandia National Laboratories
 LAMMPS development team: developers@lammps.org

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#include "atom_vec_sph.h"

#include "atom.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecSPH::AtomVecSPH(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = Atom::ATOMIC;
  mass_type = PER_TYPE;
  forceclearflag = 1;

  atom->esph_flag = 1;
  atom->rho_flag = 1;
  atom->cv_flag = 1;
  atom->vest_flag = 1;
  atom->gradv_flag = 1;
  atom->gammadot_flag = 1;
  atom->fintx_flag = 1;
  atom->finty_flag = 1;
  
  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = {"rho", "drho", "esph", "desph", "cv", "vest", "gradv", "gammadot", "fintx", "finty"};
  fields_copy = {"rho", "drho", "esph", "desph", "cv", "vest", "gradv", "gammadot", "fintx", "finty"};
  fields_comm = {"rho", "esph", "vest", "gradv", "gammadot", "fintx", "finty"};
  fields_comm_vel = {"rho", "esph", "vest", "gradv", "gammadot", "fintx", "finty"};
  fields_reverse = {"drho", "desph", "gradv", "gammadot", "fintx", "finty"};
  fields_border = {"rho", "esph", "cv", "vest", "gradv", "gammadot",  "fintx", "finty"};
  fields_border_vel = {"rho", "esph", "cv", "vest", "gradv", "gammadot", "fintx", "finty"};
  fields_exchange = {"rho", "esph", "cv", "vest", "gradv", "gammadot", "fintx", "finty"};
  fields_restart = {"rho", "esph", "cv", "vest", "gradv", "gammadot", "fintx", "finty"};
  fields_create = {"rho", "esph", "cv", "vest", "desph", "drho", "gradv", "gammadot", "fintx", "finty"};
  fields_data_atom = {"id", "type", "rho", "esph", "cv", "x"};
  fields_data_vel = {"id", "v"};

  setup_fields();
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecSPH::grow_pointers()
{
  rho = atom->rho;
  drho = atom->drho;
  esph = atom->esph;
  desph = atom->desph;
  cv = atom->cv;
  vest = atom->vest;
  gradv = atom->gradv;
  gammadot = atom->gammadot;    
  fintx = atom->fintx;
  finty = atom->finty;    
}

/* ----------------------------------------------------------------------
   clear extra forces starting at atom N
   nbytes = # of bytes to clear for a per-atom vector
------------------------------------------------------------------------- */

void AtomVecSPH::force_clear(int n, size_t nbytes)
{
  memset(&desph[n], 0, nbytes);
  memset(&drho[n], 0, nbytes);
  
  memset(&gradv[n][0],0,nbytes);
  memset(&gradv[n][1],0,nbytes);
  memset(&gradv[n][2],0,nbytes);
  memset(&gradv[n][3],0,nbytes);
  memset(&gradv[n][4],0,nbytes);
  memset(&gradv[n][5],0,nbytes);
  memset(&gradv[n][6],0,nbytes);
  memset(&gradv[n][7],0,nbytes);
  memset(&gradv[n][8],0,nbytes); 

  //memset(&gammadot[n], 0, nbytes); 
  memset(&fintx[n], 0, nbytes); 
  memset(&finty[n], 0, nbytes);     
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecSPH::create_atom_post(int ilocal)
{
  cv[ilocal] = 1.0;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecSPH::data_atom_post(int ilocal)
{
  vest[ilocal][0] = 0.0;
  vest[ilocal][1] = 0.0;
  vest[ilocal][2] = 0.0;
  desph[ilocal] = 0.0;
  drho[ilocal] = 0.0;
  
  gradv[ilocal][0] = 0.0;
  gradv[ilocal][1] = 0.0;
  gradv[ilocal][2] = 0.0;
  gradv[ilocal][3] = 0.0;
  gradv[ilocal][4] = 0.0;
  gradv[ilocal][5] = 0.0;
  gradv[ilocal][6] = 0.0;
  gradv[ilocal][7] = 0.0;
  gradv[ilocal][8] = 0.0;

  gammadot[ilocal] = 0.0;
  fintx[ilocal] = 0.0;    
  finty[ilocal] = 0.0;    
}

/* ----------------------------------------------------------------------
   assign an index to named atom property and return index
   return -1 if name is unknown to this atom style
------------------------------------------------------------------------- */

int AtomVecSPH::property_atom(const std::string &name)
{
    if (name == "rho") return 0;
    if (name == "drho") return 1;
    if (name == "esph") return 2;
    if (name == "desph") return 3;
    if (name == "cv") return 4;
    
    if (name == "gdvxx") return 5;
    if (name == "gdvyy") return 6;
    if (name == "gdvzz") return 7;
    if (name == "gdvxy") return 8;
    if (name == "gdvxz") return 9;
    if (name == "gdvyz") return 10;
    if (name == "gdvyx") return 11;
    if (name == "gdvzx") return 12;
    if (name == "gdvzy") return 13;

    if (name == "gammadot") return 14;
    if (name == "fintx") return 15;
    if (name == "finty") return 16;

/*  
  if (strcmp(name,"rho") == 0) return 0; 
  if (strcmp(name,"drho") == 0) return 1; 
  if (strcmp(name,"esph") == 0) return 2; 
  if (strcmp(name,"desph") == 0) return 3; 
  if (strcmp(name,"cv") == 0) return 4; 
  
  if (strcmp(name,"gdvxx") == 0) return 5; 
  if (strcmp(name,"gdvyy") == 0) return 6; 
  if (strcmp(name,"gdvzz") == 0) return 7; 
  if (strcmp(name,"gdvxy") == 0) return 8; 
  if (strcmp(name,"gdvxz") == 0) return 9; 
  if (strcmp(name,"gdvyz") == 0) return 10; 
  if (strcmp(name,"gdvyx") == 0) return 11; 
  if (strcmp(name,"gdvzx") == 0) return 12; 
  if (strcmp(name,"gdvzy") == 0) return 13; 
*/


  return -1;
}

/* ----------------------------------------------------------------------
   pack per-atom data into buf for ComputePropertyAtom
   index maps to data specific to this atom style
------------------------------------------------------------------------- */

void AtomVecSPH::pack_property_atom(int index, double *buf, int nvalues, int groupbit)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int n = 0;

  if (index == 0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = rho[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 1) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = drho[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 2) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = esph[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 3) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = desph[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 4) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = cv[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 5) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][0];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 6) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][1];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 7) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][2];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 8) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][3];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 9) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][4];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 10) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][5];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 11) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][6];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 12) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][7];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 13) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = gradv[i][8];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 14) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = gammadot[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
    } else if (index == 15) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = fintx[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 16) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = finty[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  }
}
