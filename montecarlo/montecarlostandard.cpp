#include "montecarlostandard.h"
#include "../matrix.h"
#include "../random.h"
#include "../utils.h"

#include <math.h>
#include <iostream>
#include <iomanip>

MonteCarloStandard::MonteCarloStandard(Config *config) :
    MonteCarlo(config),
    rank(config->rank()),
    stepLength(config->stepLength()),
    firstSample(true)
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
    m_moves = new vec2*[1];
    m_moves[0] = new vec2[nParticles];
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
    int nthMove = 0;
    if(recordMoves) {
        nthMove = nCycles * nParticles / (nMoves);
    }
    // initialisations of variational parameters and energies
    m_energy = m_energySquared = 0; localEnergy=0;

    //  initial trial position
    for (int i = 0; i < nParticles; i++) {
        for (int j=0; j < nDimensions; j++) {
            rOld[i][j] = stepLength*(ran2(idumMC)-0.5);
        }
        rNew[i] = rOld[i];
    }
    wave->initialize(rOld);
    // TODO Optimize step length by Newton's method
    // loop over monte carlo cycles
    int move = 0;
    for (cycle = 0; cycle < nCycles; cycle++){
        // new position
        for (int i = 0; i < nParticles; i++) {
            for (int j=0; j < nDimensions; j++) {
                rNew[i][j] = rOld[i][j]+stepLength*(ran2(idumMC)-0.5);
            }
            double ratio = wave->ratio(rNew[i], i);
            //            std::cout << "Ratio calculated " << std::endl;
            //            for(int i = 0; i < nParticles; i++) {
            //                std::cout << "rNew[" << i << "] = " << rNew[i] << std::endl;
            //            }
            // The Metropolis test is performed by moving one particle at the time
            if(ran2(idumMC) <= (ratio*ratio)) {
                rOld[i] = rNew[i];
                wave->acceptMove(i);
                //                std::cout << "Accepted" << std::endl;
            } else {
                rNew[i] = rOld[i]; // Move the particle back
                wave->rejectMove();
                //                std::cout << "Denied" << std::endl;
            }
            //            std::cout << "Move decided" << std::endl;
            //            for(int i = 0; i < nParticles; i++) {
            //                std::cout << "rNew[" << i << "] = " << rNew[i] << std::endl;
            //            }
            // compute local energy
            localEnergy = hamiltonian->energy(wave, rNew);
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
                if(recordMoves) {
                    if(!(cycle % nthMove)) {
                        //                    std::cout << "Recording move " << move << " @ " << cycle << std::endl;
                        m_moves[move][i] = rOld[i];
                        if(i == nParticles - 1) {
                            move++;
                        }
                    }
                }
            } else {
                checkTerminalization(localEnergy);
            }


        }   // end of loop over MC trials
    }  //  end of loop over particles
    m_energy /= (nCycles * nParticles);
    m_energySquared /= (nCycles * nParticles);
    //    std::cout << "Done sampling. Had " << terminalizationTrials << " terminalization trials with the last diff at " << diffAverage << std::endl;
}

