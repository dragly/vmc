#include "metropolishastingsmontecarlo.h"
#include <math.h>
#include "../matrix.h"
#include "../random.h"
#include "../utils.h"
#include "../config.h"

#include <iomanip>

using std::setprecision;

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
    ofstream positionFile;
    if(recordMoves) {
        positionFile.open("vmc-positions-init.dat");
    }
    wave->initialize(rOld);

    //    wave->gradient(rOld, 0, waveGradientOld);
    wave->gradient(rOld, quantumForceOld);
    quantumForceOld *= 2;
    int acceptances = 0;
    int rejections = 0;
    // loop over monte carlo cycles
    for (cycle = 0; cycle <= nCycles; cycle++){
        // new trial position
        for (int i = 0; i < nParticles; i++) {
//            wave->gradient(rOld, quantumForceOld);
//            quantumForceOld *= 2;
            for (int k=0; k < nDimensions; k++) {
                int qfIndex = i * nDimensions + k;
                rNew[i][k]= rOld[i][k] + stepLength * diffConstant * quantumForceOld[qfIndex] + 2 * diffConstant * sqrt(stepLength) * gaussianDeviate(idumMC);
            }
            wave->prepareGradient(rNew[i], i);
            wave->gradient(rNew, quantumForceNew);
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
            double ratio = wave->ratio(rNew[i], i);
            double greensRatio = exp(argSum);
            double weight = ratio*ratio * greensRatio;
            //            std::cout << ratio << std::endl;
            if(ran3(idumMC) < weight) {
                rOld[i] = rNew[i];
                wave->acceptMove(i);
                quantumForceOld = quantumForceNew;
                if(terminalized) {
                    acceptances++;
                }
                //                std::cout << "Accepted" << std::endl;
            } else {
                rNew[i] = rOld[i]; // Move the particle back
                wave->rejectMove();
                quantumForceNew = quantumForceOld;
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
                    if(!(cycle % 10)) {
                        positionFile << rOld[i][0] << " " << rOld[i][1] << std::endl;
                    }
                    if(!(cycle % nthMove)) {
                        //                    std::cout << "Recording move " << move << " @ " << cycle << std::endl;
                        m_moves[move][i] = rOld[i];
                        if(i == nParticles - 1) {
                            move++;
                        }
                    }
                }
                if(i == 0 && !(cycle % 10000) && cycle > 0) {
                    std::cout << "Cycle " << cycle << ". Current average energy is " << setprecision(16) << m_energy / (cycle * nParticles) << std::endl;
                    std::cout << "Acceptance ratio: " << (double)acceptances / (double)(rejections + acceptances) << std::endl;
                }
            } else {
                checkTerminalization(localEnergy);
            }
        }  //  end of loop over particles
    }
    if(recordMoves) {
        positionFile.close();
    }
    m_energy /= (nCycles * nParticles);
    m_energySquared /= (nCycles * nParticles);
}
