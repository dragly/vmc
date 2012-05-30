#include "diffusionwalker.h"

#include "../random.h"
#include "../wavefunction/waveslater.h"
#include "../slater/slater.h"
#include "../jastrow/jastrow.h"

DiffusionWalker::DiffusionWalker(Config *config, DiffusionWalker **otherWalkers_, int nOtherWalkers_) :
    Walker(config),
    otherWalkers(otherWalkers_),
    nWalkersMax(nOtherWalkers_),
    diffConstant(config->diffusionConstant()),
    timeStep(0.01),
    m_aliveNew(false),
    m_aliveOld(false)
{
}

void DiffusionWalker::advance(double trialEnergy) {
    Walker::advance();
    m_changeInEnergySamples = 0;
    m_energy = 0;
    for(int i = 0; i < nParticles; i++) {
        wave->gradient(rOld, quantumForceOld);
        quantumForceOld *= 2;
        // Propose move (with quantum force)
        for(int k = 0; k < nDimensions; k++) {
            // TODO per cartesian component tau?
            int qfIndex = i * nDimensions + k;
            rNew[i][k] = rOld[i][k] + timeStep * diffConstant * quantumForceOld[qfIndex] + 2 * diffConstant * sqrt(timeStep) * simpleGaussRandom(idum);
//            rNew[i][k] = rOld[i][k] + timeStep * diffConstant * quantumForceOld[qfIndex] + 2 * diffConstant * sqrt(timeStep) * randomVec[0];
//                    std::cout << "i ndim k " << i * nDimensions + k << " " << quantumForceNew->n_elem << std::endl;
//                    std::cout << "Quantum force " << j << " " << i << " " << k << " " << quantumForceNew[j][i * nDimensions + k] << std::endl;
//                    std::cout << "Old/New: " << rOld[j][i][k] << " " << rNew[j][i][k] << std::endl;
        }
        double ratio = wave->ratio(rNew[i], i);
        // Apply fixed node approximation (keep sign or reject move)
        if(ratio > 0) {
            // Compute weight function
            wave->gradient(rNew, quantumForceNew);
            quantumForceNew *= 2;
            double argSum = 0;
            for(int j = 0; j < nDimensions; j++) { // TODO figure out if it is necessary to do this with all particles (is it not zero for non-moved particles?).
                for(int k = 0; k < nDimensions; k++) {
                    int qfIndex = j * nDimensions + k;
                    double quantumForceSum = quantumForceOld[qfIndex] + quantumForceNew[qfIndex];
                    double qfPositionDiff = 0.5 * timeStep * diffConstant * (quantumForceOld[qfIndex] - quantumForceNew[qfIndex]) - (rNew[j][k] - rOld[j][k]);
                    argSum += 0.5 * quantumForceSum * qfPositionDiff;
                }
            }
            double greensRatio = exp(argSum);

            // More expensive method to calculate the Greens factor:
//            double greensRatio = 1;
//            for(int j = 0; j < nParticles; j++) {
//                for(int k = 0; k < nDimensions; k++) {
//                    int qfIndex = j * nDimensions + k;
//                    double squareArgNew = (rOld[j][k] - rNew[j][k] - tau * quantumForceNew[qfIndex] / 2);
//                    double newGreen = exp(-squareArgNew*squareArgNew / (2*tau));
//                    double squareArgOld = (rNew[j][k] - rOld[j][k] - tau * quantumForceOld[qfIndex] / 2);
//                    double oldGreen = exp(-squareArgOld*squareArgOld / (2*tau));
//                    greensRatio *= newGreen / oldGreen;
//                }
//            }

            double weight = ratio*ratio * greensRatio;
            // Accept move according to Metropolis probability
            if(weight > ran3(idum)) {
                rOld[i] = rNew[i];
                wave->acceptMove(i);
                quantumForceOld = quantumForceNew;
            } else {
                wave->rejectMove();
                rNew[i] = rOld[i];
            }
        } else {
            wave->rejectMove();
            rNew[i] = rOld[i];
        }
        // Compute branching factor PB
        localEnergyNew = hamiltonian->energy(wave, rNew);
//        localEnergyOld = hamiltonian->energy(wave, rOld);
//        double branchingFactor = exp(- timeStep * diffConstant * (0.5 * (localEnergyNew + localEnergyOld) - trialEnergy));
        double branchingFactor = exp(- timeStep * (0.5 * (localEnergyNew + localEnergyOld) - trialEnergy));
//        std::cout << aliveNew() << aliveOld() << " " << branchingFactor << " " << trialEnergy << " " << localEnergyNew << " " << localEnergyOld << std::endl;
        // Make int(PB + u) copies
        int reproductions = int(branchingFactor + ran3(idum));
        // Accumulate the energy and any observables weighted by PB
        m_energy += localEnergyNew * branchingFactor;
        m_changeInEnergySamples++;
        if(reproductions == 0 || localEnergyNew < trialEnergy - 1. / sqrt(timeStep) || localEnergyNew > trialEnergy + 1 / sqrt(timeStep)) {
            m_aliveNew = false;
        } else {
            if(reproductions > 1) {
                for(int repro = 0; repro < reproductions; repro++) {
                    for(int walkerIndex = 0; walkerIndex < nWalkersMax; walkerIndex++) {
                        DiffusionWalker *otherWalker = otherWalkers[walkerIndex];
                        // find a dead walker to ressurect
                        if(!otherWalker->aliveNew()) {
                            otherWalker->setAliveNew(true);
                            otherWalker->copyOtherWalker(this);
                            break;
                        }
                    }
                }
            }
        } // END if branching
//                std::cout << "Alive walkers: " << nWalkersAlive << std::endl;
//            std::cout << "Should make " << (int) (branchingFactor + ran3(idum)) << " copies" << std::endl;
        localEnergyOld = localEnergyNew;
    } // END for every particle
}

void DiffusionWalker::progressToNextStep()
{
    Walker::progressToNextStep();
    m_aliveOld = m_aliveNew;
}
