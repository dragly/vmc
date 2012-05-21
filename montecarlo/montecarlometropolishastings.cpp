#include "montecarlometropolishastings.h"
#include <math.h>
#include "../matrix.h"
#include "../random.h"
#include "../utils.h"
#include "../config.h"

MonteCarloMetropolisHastings::MonteCarloMetropolisHastings(Config *config) :
    MonteCarlo(config),
    rank(config->rank()),
    stepLength(config->stepLength()),
    hamiltonian(config->hamiltonian()),
    firstSample(true)
{
    // allocate matrices which contain the position of the particles
    rOld = new vec2[ nParticles];
    rNew = new vec2[ nParticles];

    quantumForceNew = zeros<vec>(nParticles * nDimensions);
    quantumForceOld = zeros<vec>(nParticles * nDimensions);
}

MonteCarloMetropolisHastings::~MonteCarloMetropolisHastings()
{
}

void MonteCarloMetropolisHastings::sample(int nCycles)
{
    m_energy = 0;
    m_energySquared = 0;
    terminalizationSum = 0;
    terminalizationNum = 1;
    double localEnergy = 0;
    double diffConstant = 0.5;
    int nthMove = 0;
    if(storeEnergies) {
        m_allEnergies = new double[nCycles];
    }
    if(recordMoves) {
        nthMove = nCycles * nParticles / (nMoves);
    }
    //  initial trial position, note calling with alpha
    for (int i = 0; i < nParticles; i++) {
        for (int j=0; j < nDimensions; j++) {
            rOld[i][j] = stepLength*(ran2(idumMC)-0.5);
        }
        rNew[i] = rOld[i];
    }
    wave->initialize(rOld);
    //    wave->gradient(rOld, 0, waveGradientOld); // TODO add particle number
    wave->gradient(rOld, 0, quantumForceOld); // TODO add particle number
    quantumForceOld *= 2;
    int acceptances = 0;
    int rejections = 0;
    stepLength = 0.01;
    // loop over monte carlo cycles
    for (int cycle = 0; cycle <= nCycles; cycle++){
        // new trial position
        for (int i = 0; i < nParticles; i++) {
            wave->gradient(rOld, 0, quantumForceOld); // TODO add particle number
            quantumForceOld *= 2;
            for (int k=0; k < nDimensions; k++) {
                int qfIndex = i * nDimensions + k;
                rNew[i][k]= rOld[i][k] +  diffConstant * quantumForceOld[qfIndex] * stepLength + sqrt(2 * diffConstant * stepLength) * simpleGaussRandom(idumMC);
            }

            // The Metropolis test is performed by moving one particle at the time
            wave->gradient(rNew, 0, quantumForceNew); // TODO add particle number
            quantumForceNew *= 2;
            double argSum = 0;
            for(int j = 0; j < nParticles; j++) {
                for(int k = 0; k < nDimensions; k++) {
                    int qfIndex = j * nDimensions + k;
                    double quantumForceSum = quantumForceOld[qfIndex] + quantumForceNew[qfIndex];
                    double qfPositionDiff = 0.5 * diffConstant * stepLength * (quantumForceOld[qfIndex] - quantumForceNew[qfIndex]) - (rNew[j][k] - rOld[j][k]);
                    argSum += 0.5 * quantumForceSum * qfPositionDiff;
                }
            }
            double greensRatio = exp(argSum);

            double ratio = wave->ratio(rNew[i], i);
            double weight = ratio*ratio * greensRatio;
            if(ran2(idumMC) <= weight) {
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
        }  //  end of loop over particles
    }
//    std::cout << "Acceptance ratio: " << (double)acceptances / (double)(rejections + acceptances) << std::endl;
    m_energy /= (nCycles * nParticles);
    m_energySquared /= (nCycles * nParticles);
}
