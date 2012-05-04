#include <math.h>

#include "montecarlostandard.h"
#include "../matrix.h"
#include "../random.h"
#include "../utils.h"

MonteCarloStandard::MonteCarloStandard(Config *config) :
    MonteCarlo(config),
    rank(config->rank()),
    step_length(config->stepLength()),
    wave(config->wave())
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
    wave->init(rOld);
    // TODO Optimize step length by Newton's method
    // loop over monte carlo cycles
    for (cycle = 0; cycle < nCycles; cycle++){
        // new position
        for (int i = 0; i < nParticles; i++) {
            for (int j=0; j < nDimensions; j++) {
                rNew[i][j] = rOld[i][j]+step_length*(ran2(idum)-0.5);
            }
            double ratio = wave->ratio(rNew);
            // The Metropolis test is performed by moving one particle at the time
            if(ran2(idum) <= (ratio*ratio)) {
                rOld[i] = rNew[i];
                wave->acceptEvaluation();
            } else {
                rNew[i] = rOld[i]; // Move the particle back
            }
            // compute local energy
            localEnergy = config->hamiltonian()->energy(wave, rOld);
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

