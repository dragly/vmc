#include <math.h>

#include "montecarlostandard.h"
#include "matrix.h"
#include "random.h"
#include "utils.h"

MonteCarloStandard::MonteCarloStandard(WaveFunction *wave, Hamiltonian* hamiltonian, int number_particles, int dimension, double charge, int rank, double step_length) :
    MonteCarlo(wave, hamiltonian),
    number_particles(number_particles),
    dimension(dimension),
    charge(charge),
    rank(rank),
    step_length(step_length)
{
}

void MonteCarloStandard::sample(int max_variations, int number_cycles, double *cumulative_e, double *cumulative_e2, double *all_energies)
{
    int cycles, variate, i, j, k;
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
        energy = energy2 = 0; delta_e=0;
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
                // TODO: Optimize this by removing the if-test. Profile first!
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
            delta_e = hamiltonian->energy(wave, r_old, alpha, wfold);
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
        alpha += 0.1;
    }    // end of loop over variational  steps
    free_matrix((void **) r_old); // free memory
    free_matrix((void **) r_new); // free memory
}
