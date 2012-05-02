#include <math.h>

#include "montecarlostandard.h"
#include "../matrix.h"
#include "../random.h"
#include "../utils.h"

MonteCarloStandard::MonteCarloStandard(Config *config) :
    MonteCarlo(config),
    rank(config->rank()),
    step_length(config->stepLength())
{
    // allocate matrices which contain the position of the particles
    rOld = new vec2[ nParticles];
    rNew = new vec2[ nParticles];
    for (int i = 0; i < nParticles; i++) {
        for (int j=0; j < nDimensions; j++) {
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

void MonteCarloStandard::sample(int nCycles)
{
    double wfnew = 0;
    double wfold = 0;
    m_energy = 0;
    m_energySquared = 0;
    terminalizationSum = 0;
    terminalizationNum = 1;
    prevTerminalizationAverage = 99999999;
    double localEnergy = 0;
    if(storeEnergies) {
        m_allEnergies = new double[nCycles];
    }
    // initialisations of variational parameters and energies
    m_energy = m_energySquared = 0; localEnergy=0;
    //  initial trial position, note calling with alpha
    for (int i = 0; i < nParticles; i++) {
        for (int j=0; j < nDimensions; j++) {
            rOld[i][j] = step_length*(ran2(idum)-0.5);
        }
    }
    wfold = config->wave()->wave(rOld);
    // loop over monte carlo cycles
    for (cycle = 0; cycle < nCycles; cycle++){
        // new position
        for (int i = 0; i < nParticles; i++) {
            for (int j=0; j < nDimensions; j++) {
                rNew[i][j] = rOld[i][j]+step_length*(ran2(idum)-0.5);
            }
            // TODO Optimize MonteCarloStandard by removing the if-test. Profile first!
            //  for the other particles we need to set the position to the old position since
            //  we move only one particle at the time
            for (int k = 0; k < nParticles; k++) {
                if ( k != i) {
                    for (int l=0; l < nDimensions; l++) {
                        rNew[k][l] = rOld[k][l];
                    }
                }
            }
            wfnew = config->wave()->wave(rNew);
            // The Metropolis test is performed by moving one particle at the time
            if(ran2(idum) <= wfnew*wfnew/(wfold*wfold)) {
                for (int l=0; l < nDimensions; l++) {
                    rOld[i][l]=rNew[i][l];
                }
                wfold = wfnew;
            }
            // compute local energy
            localEnergy = config->hamiltonian()->energy(config->wave(), rOld);
            // save all energies on last variate
            //        if(variate==max_variations){
            if(terminalized) {
                if(storeEnergies) {
                    m_allEnergies[cycle] = localEnergy;
                }
                //        }
                // update energies
                m_energy += localEnergy;
                m_energySquared += localEnergy*localEnergy;
            } else {
                checkTerminalization(localEnergy);
            }
        }   // end of loop over MC trials
    }  //  end of loop over particles
    m_energy /= (nCycles * nParticles);
    m_energySquared /= (nCycles * nParticles);
    std::cout << "Done sampling. Had " << terminalizationTrials << " terminalization trials with the last diff at " << diffAverage << std::endl;
}

