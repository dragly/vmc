#include "metropolishastingsmontecarlo.h"
#include <math.h>
#include "../matrix.h"
#include "../random.h"
#include "../utils.h"
#include "../config.h"

MetropolisHastingsMonteCarlo::MetropolisHastingsMonteCarlo(Config *config) :
    MonteCarlo(config),
    myRank(config->myRank()),
    hamiltonian(config->hamiltonian()),
    firstSample(true),
    diffConstant(config->diffusionConstant())
{
    quantumForceNew = zeros<vec>(nParticles * nDimensions);
    quantumForceOld = zeros<vec>(nParticles * nDimensions);
}

MetropolisHastingsMonteCarlo::~MetropolisHastingsMonteCarlo()
{
}

void MetropolisHastingsMonteCarlo::sample(int nCycles)
{
    m_energy = 0;
    m_energySquared = 0;
    terminalizationSum = 0;
    terminalizationNum = 1;
    double localEnergy = 0;
    int nthMove = 0;
    if(storeEnergies) {
        m_allEnergies = new double[nCycles];
    }
    if(recordMoves) {
        nthMove = nCycles * nParticles / (nMoves);
    }
    //  initial trial position, note calling with alpha
//    for (int i = 0; i < nParticles; i++) {
//        for (int j=0; j < nDimensions; j++) {
//            rOld[i][j] = stepLength*(ran3(idumMC)-0.5);
//        }
//        rNew[i] = rOld[i];
//    }
    wave->initialize(rOld);
    //    wave->gradient(rOld, 0, waveGradientOld);
    wave->gradient(rOld, quantumForceOld);
    quantumForceOld *= 2;
    m_variationalGradient = zeros(2); // TODO make parameter-specific
    int acceptances = 0;
    int rejections = 0;
    // loop over monte carlo cycles
    for (cycle = 0; cycle <= nCycles; cycle++){
        // new trial position
        for (int i = 0; i < nParticles; i++) {
            wave->gradient(rOld, quantumForceOld);
            quantumForceOld *= 2;
            for (int k=0; k < nDimensions; k++) {
                int qfIndex = i * nDimensions + k;
                rNew[i][k]= rOld[i][k] +  diffConstant * quantumForceOld[qfIndex] * stepLength + 2 * diffConstant * sqrt(stepLength) * simpleGaussRandom(idumMC);
            }

            // The Metropolis test is performed by moving one particle at the time
            wave->gradient(rNew, quantumForceNew);
            quantumForceNew *= 2;
            double argSum = 0;
            for(int k = 0; k < nDimensions; k++) {
                int qfIndex = i * nDimensions + k;
                double quantumForceSum = quantumForceOld[qfIndex] + quantumForceNew[qfIndex];
                double qfPositionDiff = 0.5 * diffConstant * stepLength * (quantumForceOld[qfIndex] - quantumForceNew[qfIndex]) - (rNew[i][k] - rOld[i][k]);
                argSum += 0.5 * quantumForceSum * qfPositionDiff;
            }
            double greensRatio = exp(argSum);

            double ratio = wave->ratio(rNew[i], i);
            double weight = ratio*ratio * greensRatio;
//            std::cout << ratio << std::endl;
            if(ran3(idumMC) < weight) {
                rOld[i] = rNew[i];
                wave->acceptMove(i);
                if(terminalized) {
                    acceptances++;
                }
                //                std::cout << "Accepted" << std::endl;
            } else {
                rNew[i] = rOld[i]; // Move the particle back
                wave->rejectMove();
                if(terminalized) {
                    rejections++;
                }
            }
            localEnergy = hamiltonian->energy(wave, rOld);
            if(sampleVariationalGradient) {
                m_variationalGradient = m_variationalGradient + wave->variationalGradient();
            }
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
                if(i == 0 && !(cycle % 10000) && cycle > 0) {
                    std::cout << "Energy from last 10000 cycles was " << m_energy / (cycle * nParticles) << std::endl;
                }
            } else {
                checkTerminalization(localEnergy);
            }
        }  //  end of loop over particles
    }
    std::cout << "Acceptance ratio: " << (double)acceptances / (double)(rejections + acceptances) << std::endl;
    m_energy /= (nCycles * nParticles);
    m_energySquared /= (nCycles * nParticles);
    m_variationalGradient = m_variationalGradient / (nCycles * nParticles);
}
