#include "standardmontecarlo.h"
#include "../matrix.h"
#include "../random.h"
#include "../utils.h"

#include <math.h>
#include <iostream>
#include <iomanip>

using std::setprecision;

StandardMonteCarlo::StandardMonteCarlo(Config *config) :
    MonteCarlo(config),
    myRank(config->myRank()),
    firstSample(true)
{
    m_allEnergies = new double[1]; // dummy array that can be deleted when sampling with energy storage
    m_moves = new vec2*[1];
    m_moves[0] = new vec2[nParticles];
}

StandardMonteCarlo::~StandardMonteCarlo()
{
}
void StandardMonteCarlo::sample(int nSamples)
{
    m_energy = 0;
    m_energySquared = 0;
    terminalizationSum = 0;
    terminalizationNum = 1;
    prevTerminalizationAverage = 99999999;
    double localEnergy = 0;
    if(storeEnergies) {
        // TODO delete m_allEnergies first
        m_allEnergies = new double[nSamples];
    }
    int nthMove = 0;
    if(recordMoves) {
        nthMove = nSamples * nParticles / (nMoves);
    }
    // initialisations of variational parameters and energies
    m_energy = m_energySquared = 0; localEnergy=0;

    //  initial trial position
//    if(firstSample) {
//        for (int i = 0; i < nParticles; i++) {
//            for (int j=0; j < nDimensions; j++) {
//                rOld[i][j] = stepLength*(ran3(idumMC)-0.5);
//            }
//            rNew[i] = rOld[i];
//        }
////        firstSample = false;
//    }
    wave->initialize(rOld);
    int acceptances = 0;
    int rejections = 0;
    // TODO Optimize step length by Newton's method
    // loop over monte carlo cycles
    for (cycle = 0; cycle < nSamples; cycle++){
        // new position
        for (int i = 0; i < nParticles; i++) {
            for (int j=0; j < nDimensions; j++) {
                rNew[i][j] = rOld[i][j] + stepLength*(ran3(idumMC)-0.5);
            }
            wave->prepareGradient(rNew[i],i);
            double ratio = wave->ratio(rNew[i], i);
            //            std::cout << "Ratio calculated " << std::endl;
            //            for(int i = 0; i < nParticles; i++) {
            //                std::cout << "rNew[" << i << "] = " << rNew[i] << std::endl;
            //            }
            // The Metropolis test is performed by moving one particle at the time
            if(ran3(idumMC) <= (ratio*ratio)) {
                rOld[i] = rNew[i];
                wave->acceptMove(i);
                acceptances++;
                //                std::cout << "Accepted" << std::endl;
            } else {
                rNew[i] = rOld[i]; // Move the particle back
                wave->rejectMove();
                rejections++;
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
                    if(!(cycle % nSamples / int(1e9 / nParticles))) { // store 10^x positions
                        writePositionToFile(rOld[i]);
                    }
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
                if(terminalized) {
                    i = -1;
                }
            }

            if(outputEnergies && i == nParticles - 1 && !(cycle % 10000) && cycle > 0) {
                std::cout << "Cycle " << cycle << ". Current average energy is " << setprecision(16) << m_energy / ((cycle + 1)* nParticles) << ". Acceptance ratio: " << setprecision(6) << (double)acceptances / (double)(rejections + acceptances) << " local energy " << setprecision(10) << localEnergy << std::endl;
                hamiltonian->outputTotals();
            }

        }   // end of loop over MC trials
    }  //  end of loop over particles
    m_energy /= ((nSamples) * nParticles);
    m_energySquared /= ((nSamples) * nParticles);
    //    std::cout << "Done sampling. Had " << terminalizationTrials << " terminalization trials with the last diff at " << diffAverage << std::endl;
}

