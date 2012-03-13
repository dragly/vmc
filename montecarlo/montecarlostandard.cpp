#include <math.h>

#include "montecarlostandard.h"
#include "../matrix.h"
#include "../random.h"
#include "../utils.h"

MonteCarloStandard::MonteCarloStandard(Config *config) :
    MonteCarlo(config),
    m_nParticles(config->nParticles()),
    m_nDimensions(config->nDimensions()),
    charge(config->charge()),
    rank(config->rank()),
    step_length(config->stepLength())
{
    // every node has its own seed for the random numbers
    idum = -1-rank;
    // allocate matrices which contain the position of the particles
    r_old = (double **) matrix( m_nParticles, m_nDimensions, sizeof(double));
    r_new = (double **) matrix( m_nParticles, m_nDimensions, sizeof(double));
    for (int i = 0; i < m_nParticles; i++) {
        for (int j=0; j < m_nDimensions; j++) {
            r_old[i][j] = r_new[i][j] = 0;
        }
    }
}

MonteCarloStandard::~MonteCarloStandard()
{
    free_matrix((void **) r_old); // free memory
    free_matrix((void **) r_new); // free memory
}

void MonteCarloStandard::sample(int nCycles, double *energies, double *allEnergies)
{
    double wfnew = 0;
    double wfold = 0;
    double energy = 0;
    double energy2 = 0;
    double delta_e = 0;
    // initialisations of variational parameters and energies
    energy = energy2 = 0; delta_e=0;
    //  initial trial position, note calling with alpha
    for (int i = 0; i < m_nParticles; i++) {
        for (int j=0; j < m_nDimensions; j++) {
            r_old[i][j] = step_length*(ran2(&idum)-0.5);
        }
    }
    wfold = m_config->wave()->wave(r_old);
    // loop over monte carlo cycles
    for (int cycle = 1; cycle <= nCycles; cycle++){
        // new position
        for (int i = 0; i < m_nParticles; i++) {
            for (int j=0; j < m_nDimensions; j++) {
                r_new[i][j] = r_old[i][j]+step_length*(ran2(&idum)-0.5);
            }
            // TODO Optimize MonteCarloStandard by removing the if-test. Profile first!
            //  for the other particles we need to set the position to the old position since
            //  we move only one particle at the time
            for (int k = 0; k < m_nParticles; k++) {
                if ( k != i) {
                    for (int l=0; l < m_nDimensions; l++) {
                        r_new[k][l] = r_old[k][l];
                    }
                }
            }
            wfnew = m_config->wave()->wave(r_new);
            // The Metropolis test is performed by moving one particle at the time
            if(ran2(&idum) <= wfnew*wfnew/(wfold*wfold)) {
                for (int l=0; l < m_nDimensions; l++) {
                    r_old[i][l]=r_new[i][l];
                }
                wfold = wfnew;
            }
        }  //  end of loop over particles
        // compute local energy
        delta_e = m_config->hamiltonian()->energy(m_config->wave(), r_old);
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
