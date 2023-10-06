/******************************************************************************/
/*                                                    			              */
/*  A Molecular Dynamics simulation of Lennard-Jones particles                */
/*                                                                            */
/*	This code is part of the course "Particle-based Simulations"              */
/*  taught at Eindhoven University of Technology.                             */
/*  No part of this code may be reproduced without permission of the author:  */
/*  Dr. Ir. E.A.J.F. Peters                                                   */
/*                                                                            */
/*  Dr. Ir. J.T. Padding:    version 1.1, 30/1/2013                           */
/*  Jeroen Hofman:           version 1.2, 28/7/2015                           */
/*  Dr. Ir. E.A.J.F. Peters: version 4.0, 18/9/2018    			              */
/******************************************************************************/

/*
 * For the 2023 PBS assigment the code needs to be extended
 *
 * -Implement a Berendsen thermostat in dynamics.c
 * -Implement bonds in initialise_bonds in file initialise.c
 * -Initialize vectors.type such that particles get the proper type 
 * -Implement the bonded and non-bonded force in force.c. (Make the forces type dependent)
 * -Change the particle position initialisation such that it takes into account bond lengths and angles
 * -Implement the needed on-the-fly data analysis
 * 
 */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "setparameters.h"
#include "initialise.h"
#include "nbrlist.h"
#include "forces.h"
#include "dynamics.h"
#include "memory.h"
#include "fileoutput.h"
#include "RadDist.h"

/**
 * @brief main The main of the MD code. After initialization,
 * a velocity-Verlet scheme is executed for a specified number of time steps.
 *
 * @return int 0 if successful
 */

void compute_velocity_distribution(struct Vectors* p_vectors, struct Parameters* p_parameters, double* histogram, int num_bins, double v_max);

int main(void)
{
    struct Vectors vectors;
    struct Parameters parameters;
    struct Nbrlist nbrlist;
    size_t step;
    double Ekin, Epot, time;
    FILE *fpt2;

    set_parameters(&parameters);
    alloc_memory(&parameters, &vectors, &nbrlist);
    if (parameters.load_restart == 1)
    {
        load_restart(&parameters, &vectors);
        initialise_structure(&parameters, &vectors, &nbrlist);
        step = 0;
        time = 0.0;
    }
    else
        initialise(&parameters, &vectors, &nbrlist, &step, &time);
    build_nbrlist(&parameters, &vectors, &nbrlist);
    Epot = calculate_forces(&parameters, &nbrlist, &vectors,step,fpt2);
    record_trajectories_pdb(1, &parameters, &vectors, time);

    FILE *fpt;
    fpt = fopen("Save_Energy.csv","w+");
    fpt2 = fopen("radDist2.csv","w+");
    fprintf(fpt,"Time_step, Time, Epot, Ekin, Etot\n");

    // Add histogram related variables
    int num_bins = 100;  
    double v_max = 0.018; 

    double* velocity_histogram = (double*)malloc(num_bins * sizeof(double));
    if (velocity_histogram == NULL) {
        perror("Failed to allocate memory for velocity histogram");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < num_bins; i++) {
        velocity_histogram[i] = 0.0;
    }

    int start_collecting_histogram = 0; 

    while (step < parameters.num_dt_steps)  // start of the velocity-Verlet loop
    {
        step++;
        time += parameters.dt;
        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors);
        thermostat(&parameters, &vectors, Ekin);
        update_positions(&parameters, &nbrlist, &vectors);
        boundary_conditions(&parameters, &vectors);
        update_nbrlist(&parameters, &vectors, &nbrlist);
        Epot = calculate_forces(&parameters, &nbrlist, &vectors,step,fpt2);
        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors);

        printf("Step %zu, Time %f, Epot %f, Ekin %f, Etot %f\n", step, time, Epot, Ekin, Epot + Ekin);
        fprintf(fpt,"%zu, %f, %f, %f, %f\n", step, time, Epot, Ekin, Epot + Ekin);

        if (step > start_collecting_histogram) {
            compute_velocity_distribution(&vectors, &parameters, velocity_histogram, num_bins, v_max);
        }

        if (step % parameters.num_dt_pdb == 0) record_trajectories_pdb(0, &parameters, &vectors, time);
        if (step % parameters.num_dt_restart == 0) save_restart(&parameters,&vectors); 
    }

    fclose(fpt);

    // Normalize histogram
    double total_count = 0.0;
    for (int i = 0; i < num_bins; i++) {
        total_count += velocity_histogram[i];
    }
    for (int i = 0; i < num_bins; i++) {
        velocity_histogram[i] /= total_count;
    }

    // Output histogram to file
    FILE* hist_file = fopen("Velocity_Distribution.csv", "w+");
    if (hist_file == NULL) {
        perror("Failed to open histogram file");
        exit(EXIT_FAILURE);
    }

    fprintf(hist_file, "Velocity, Probability\n");
    double dv = v_max / num_bins;
    for (int i = 0; i < num_bins; i++) {
        fprintf(hist_file, "%f, %f\n", i * dv, velocity_histogram[i]);
    }

    fclose(hist_file);

    save_restart(&parameters, &vectors);
    free_memory(&vectors, &nbrlist);

    free(velocity_histogram);

    return 0;
}

void compute_velocity_distribution(struct Vectors* p_vectors, struct Parameters* p_parameters, double* histogram, int num_bins, double v_max) {
    for (size_t i = 0; i < p_parameters->num_part; i++) {
        double v = sqrt(p_vectors->v[i].x * p_vectors->v[i].x + p_vectors->v[i].y * p_vectors->v[i].y + p_vectors->v[i].z * p_vectors->v[i].z);
        int bin = (int)(v / v_max * num_bins);
        if (bin >= 0 && bin < num_bins) {
            histogram[bin]++;
        }
    }
}