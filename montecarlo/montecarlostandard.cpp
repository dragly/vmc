#include <math.h>

#include "montecarlostandard.h"
#include "../matrix.h"
#include "../random.h"
#include "../utils.h"

MonteCarloStandard::MonteCarloStandard(Config *config) :
    MonteCarlo(config),
    m_nParticles(config->nParticles()),
    m_nDimensions(config->nDimensions()),
    rank(config->rank()),
    step_length(config->stepLength())
{
    // allocate matrices which contain the position of the particles
    rOld = new vec2[ m_nParticles];
    rNew = new vec2[ m_nParticles];
    for (int i = 0; i < m_nParticles; i++) {
        for (int j=0; j < m_nDimensions; j++) {
            rOld[i][j] = rNew[i][j] = 0;
        }
    }
    m_allEnergies = new double[1]; // dummy array that can be deleted when sampling with energy storage
}

MonteCarloStandard::~MonteCarloStandard()
{
    delete [] rOld;
    delete [] rNew;
}

void MonteCarloStandard::sample(int nCycles, bool storeEnergies)
{
    double wfnew = 0;
    double wfold = 0;
    m_energy = 0;
    m_energySquared = 0;
    double delta_e = 0;
    if(storeEnergies) {
        m_allEnergies = new double[nCycles];
    }
    // initialisations of variational parameters and energies
    m_energy = m_energySquared = 0; delta_e=0;
    //  initial trial position, note calling with alpha
    for (int i = 0; i < m_nParticles; i++) {
        for (int j=0; j < m_nDimensions; j++) {
            rOld[i][j] = step_length*(ran2(idum)-0.5);
        }
    }
    wfold = m_config->wave()->wave(rOld);
    // loop over monte carlo cycles
    for (int cycle = 1; cycle <= nCycles; cycle++){
        // new position
        for (int i = 0; i < m_nParticles; i++) {
            for (int j=0; j < m_nDimensions; j++) {
                rNew[i][j] = rOld[i][j]+step_length*(ran2(idum)-0.5);
            }
            // TODO Optimize MonteCarloStandard by removing the if-test. Profile first!
            //  for the other particles we need to set the position to the old position since
            //  we move only one particle at the time
            for (int k = 0; k < m_nParticles; k++) {
                if ( k != i) {
                    for (int l=0; l < m_nDimensions; l++) {
                        rNew[k][l] = rOld[k][l];
                    }
                }
            }
            wfnew = m_config->wave()->wave(rNew);
            // The Metropolis test is performed by moving one particle at the time
            if(ran2(idum) <= wfnew*wfnew/(wfold*wfold)) {
                for (int l=0; l < m_nDimensions; l++) {
                    rOld[i][l]=rNew[i][l];
                }
                wfold = wfnew;
            }
        }  //  end of loop over particles
        // compute local energy
        delta_e = m_config->hamiltonian()->energy(m_config->wave(), rOld);
        // save all energies on last variate
        //        if(variate==max_variations){
        if(storeEnergies) {
            m_allEnergies[cycle] = delta_e;
        }
        //        }
        // update energies
        m_energy += delta_e;
        m_energySquared += delta_e*delta_e;
    }   // end of loop over MC trials
}
