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
    // every node has its own seed for the random numbers
    idum = -1-rank;
    // allocate matrices which contain the position of the particles
    r_old = (double **) matrix( number_particles, dimension, sizeof(double));
    r_new = (double **) matrix( number_particles, dimension, sizeof(double));
    for (int i = 0; i < number_particles; i++) {
        for (int j=0; j < dimension; j++) {
            r_old[i][j] = r_new[i][j] = 0;
        }
    }
}

void MonteCarloStandard::sample(int nCycles, double *energies, double *allEnergies)
{
    int cycle, i, j, k;
    long idum;
    double wfnew;
    double wfold;
    double energy;
    double energy2;
    double delta_e;
    // initialisations of variational parameters and energies
    energy = energy2 = 0; delta_e=0;
    //  initial trial position, note calling with alpha
    for (i = 0; i < number_particles; i++) {
        for ( j=0; j < dimension; j++) {
            r_old[i][j] = step_length*(ran2(&idum)-0.5);
        }
    }
    wfold = wave->wave(r_old);
    // loop over monte carlo cycles
    for (cycle = 1; cycle <= nCycles; cycle++){
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
            wfnew = wave->wave(r_new);
            // The Metropolis test is performed by moving one particle at the time
            if(ran2(&idum) <= wfnew*wfnew/wfold/wfold ) {
                for ( j=0; j < dimension; j++) {
                    r_old[i][j]=r_new[i][j];
                }
                wfold = wfnew;
            }
        }  //  end of loop over particles
        // compute local energy
        delta_e = hamiltonian->energy
                (wave, r_old);
        // save all energies on last variate
        //        if(variate==max_variations){
        allEnergies[cycle] = delta_e;
        //        }
        // update energies
        energy += delta_e;
        energy2 += delta_e*delta_e;
    }   // end of loop over MC trials
    // return the energy and the energy squared
    energies[0] = energy;
    energies[1] = energy2;
}

MonteCarloStandard::~MonteCarloStandard()
{
    free_matrix((void **) r_old); // free memory
    free_matrix((void **) r_new); // free memory
}
