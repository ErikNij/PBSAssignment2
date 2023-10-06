#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

void set_parameters(struct Parameters *p_parameters)
/* Set the parameters of this simulation */
{
// The parameters first 5 parameters are only used for demonstration puprposes
  p_parameters->T = 298.15;       // simulation at 25 degrees celcius
  p_parameters->kT = Kb*p_parameters->T; //2.48E-4;      //(Angstrong^2 u )/(fs^2)                                   //thermal energy
  
  p_parameters->mass = 16.04;      // u                                 //mass of a particle
  
  
  p_parameters->epsilon = 148*Kb;//1.23E-4; //(Amstrong^2 u)/fs^2; 
  p_parameters->epsilonArray[0] = 46.0*Kb;    //epsilon for CH3
  p_parameters->epsilonArray[1] = 98.0*Kb;    //epsilon for CH2
  p_parameters->epsilonArray[2] = 46.0*Kb;    //epsilon for CH2
  
  p_parameters->sigma = 3.73;      // Angstrong                                //LJ particle diameter
  p_parameters->sigmaArray[0] = 3.95;         //sigma for CH3
  p_parameters->sigmaArray[1] = 3.75;         //sigma for CH2
  p_parameters->sigmaArray[2] = 3.95;         //sigma for CH2

  p_parameters->massArray[0] = 12.0 + 3.0;
  p_parameters->massArray[1] = 12.0 + 2.0;
  p_parameters->massArray[2] = 12.0 + 3.0;


   p_parameters->r_0 = 1.54;                // Angstrom   Equilibrium diameter
   p_parameters->k_b = Kb*3.19e5;            // (u/(fs^2*K))   Angular stiffness
   p_parameters->k_t = Kb*62e3;
   p_parameters->theta_0 = 1.98968;


// The parameters below control core functionalities of the code, but many values will need to be changed
  p_parameters->num_part = 999;                            //number of particles
  p_parameters->num_dt_steps = 200;                        //number of time steps
  p_parameters->exclude_12_nb = 1;                          // 1-2 connected atoms exluded from non-bonded interactions 
  p_parameters->exclude_13_nb = 1;                          // 1-3 connected atoms exluded from non-bonded interactions    
  p_parameters->dt = 0.001;                                  //integration time step
  p_parameters->tau = 10;                                // typical 0.1 picosecond = 10 fentosecond
  //p_parameters->L = (struct Vec3D){14.938, 14.938, 14.938}; //box size                                 
  p_parameters->L = (struct Vec3D){60.0, 60.0, 60.0};       // box size for question 4 to 12  
    p_parameters->r_cut = 2.5;                              //cut-off distance used for neigbor list
  p_parameters->r_shell = 0.4;                              //shell thickness for neighbor list
  p_parameters->num_dt_pdb = 100;                           //number of time steps in between pdb outputs
  strcpy(p_parameters->filename_pdb, "trajectories");       //filename (without extension) for pdb file
  p_parameters->rescale_output = 1;                         //factor used to rescale output lengthscale (Most visualisation programs identify bonds based on distances of order 1)
  p_parameters->load_restart = 0;                           //if equal 1 restart file is loaded
  strcpy(p_parameters->restart_in_filename, "restart.dat"); //filename for loaded restart file
  p_parameters->num_dt_restart = 1000;                      // number of time steps between saves
  strcpy(p_parameters->restart_out_filename, "restart.dat");//filename for saved restart file


  if (p_parameters->r_cut > p_parameters->L.x / 2.0)
    fprintf(stderr, "Warning! r_cut > Lx/2");
  if (p_parameters->r_cut > p_parameters->L.y / 2.0)
    fprintf(stderr, "Warning! r_cut > Ly/2");
  if (p_parameters->r_cut > p_parameters->L.z / 2.0)
    fprintf(stderr, "Warning! r_cut > Lz/2");
}