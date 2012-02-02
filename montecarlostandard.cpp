#include <math.h>

#include "montecarlostandard.h"
#include "matrix.h"
#include "random.h"

// the step length and its squared inverse for the second derivative
#define h 0.001
#define h2 1000000

MonteCarloStandard::MonteCarloStandard(WaveFunction *wave, int number_particles, int dimension, double charge, int rank, double step_length) :
    MonteCarlo(wave),
    number_particles(number_particles),
    dimension(dimension),
    charge(charge),
    step_length(step_length)
{
}

void MonteCarloStandard::sample(int max_variations, int number_cycles, double *cumulative_e, double *cumulative_e2, double *all_energies)
{
    int cycles, variate, accept, dim, i, j, k;
    long idum;
    double wfnew, wfold, alpha, energy, energy2, delta_e;
    double **r_old, **r_new;
    alpha = 0.5*charge;

    // every node has its own seed for the random numbers
    idum = -1-rank;
    // allocate matrices which contain the position of the particles
    r_old = (double **) matrix( number_particles, dimension, sizeof(double));
    r_new = (double **) matrix( number_particles, dimension, sizeof(double));
    for (i = 0; i < number_particles; i++) {
        for ( j=0; j < dimension; j++) {
            r_old[i][j] = r_new[i][j] = 0;
        }
    }
    // loop over variational parameters
    for (variate=1; variate <= max_variations; variate++){
        // initialisations of variational parameters and energies
        alpha += 0.1;
        energy = energy2 = 0; accept =0; delta_e=0;
        //  initial trial position, note calling with alpha
        for (i = 0; i < number_particles; i++) {
            for ( j=0; j < dimension; j++) {
                r_old[i][j] = step_length*(ran2(&idum)-0.5);
            }
        }
        wfold = wave->wave(r_old, alpha);
        // loop over monte carlo cycles
        for (cycles = 1; cycles <= number_cycles; cycles++){
            // new position
            for (i = 0; i < number_particles; i++) {
                for ( j=0; j < dimension; j++) {
                    r_new[i][j] = r_old[i][j]+step_length*(ran2(&idum)-0.5);
                }
                //  for the other particles we need to set the position to the old position since
                //  we move only one particle at the time
                for (k = 0; k < number_particles; k++) {
                    if ( k != i) {
                        for ( j=0; j < dimension; j++) {
                            r_new[k][j] = r_old[k][j];
                        }
                    }
                }
                wfnew = wave->wave(r_new, alpha);
                // The Metropolis test is performed by moving one particle at the time
                if(ran2(&idum) <= wfnew*wfnew/wfold/wfold ) {
                    for ( j=0; j < dimension; j++) {
                        r_old[i][j]=r_new[i][j];
                    }
                    wfold = wfnew;
                }
            }  //  end of loop over particles
            // compute local energy
            delta_e = local_energy(r_old, alpha, wfold);
            // save all energies on last variate
            if(variate==max_variations){
                all_energies[cycles] = delta_e;
            }
            // update energies
            energy += delta_e;
            energy2 += delta_e*delta_e;
        }   // end of loop over MC trials
        // update the energy average and its squared
        cumulative_e[variate] = energy;
        cumulative_e2[variate] = energy2;
    }    // end of loop over variational  steps
    free_matrix((void **) r_old); // free memory
    free_matrix((void **) r_new); // free memory
}

// Function to calculate the local energy with num derivative

double MonteCarloStandard::local_energy(double **r, double alpha, double wfold)
{
    int i, j , k;
    double e_local, wfminus, wfplus, e_kinetic, e_potential, r_12,
            r_single_particle;
    double **r_plus, **r_minus;

    // allocate matrices which contain the position of the particles
    // the function matrix is defined in the progam library
    r_plus = (double **) matrix( number_particles, dimension, sizeof(double));
    r_minus = (double **) matrix( number_particles, dimension, sizeof(double));
    for (i = 0; i < number_particles; i++) {
        for ( j=0; j < dimension; j++) {
            r_plus[i][j] = r_minus[i][j] = r[i][j];
        }
    }
    // compute the kinetic energy
    e_kinetic = 0;
    for (i = 0; i < number_particles; i++) {
        for (j = 0; j < dimension; j++) {
            r_plus[i][j] = r[i][j]+h;
            r_minus[i][j] = r[i][j]-h;
            wfminus = wave->wave(r_minus, alpha);
            wfplus  = wave->wave(r_plus, alpha);
            e_kinetic -= (wfminus+wfplus-2*wfold);
            r_plus[i][j] = r[i][j];
            r_minus[i][j] = r[i][j];
        }
    }
    // include electron mass and hbar squared and divide by wave function
    e_kinetic = 0.5*h2*e_kinetic/wfold;
    // compute the potential energy
    e_potential = 0;
    // contribution from electron-proton potential
    for (i = 0; i < number_particles; i++) {
        r_single_particle = 0;
        for (j = 0; j < dimension; j++) {
            r_single_particle += r[i][j]*r[i][j];
        }
        e_potential -= charge/sqrt(r_single_particle);
    }
    // contribution from electron-electron potential
    for (i = 0; i < number_particles-1; i++) {
        for (j = i+1; j < number_particles; j++) {
            r_12 = 0;
            for (k = 0; k < dimension; k++) {
                r_12 += (r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
            }
            e_potential += 1/sqrt(r_12);
        }
    }
    free_matrix((void **) r_plus); // free memory
    free_matrix((void **) r_minus);
    e_local = e_potential+e_kinetic;
    return e_local;
}
