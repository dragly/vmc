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
    tau(config->tau()),
    m_aliveNew(false),
    m_aliveOld(false)
{
}

void DiffusionWalker::advance(double trialEnergy) {
    Walker::advance();
    m_changeInEnergySamples = 0;
    m_energy = 0;
    for(int i = 0; i < nParticles; i++) {

        // Propose move (with quantum force)
        for(int k = 0; k < nDimensions; k++) {
            // TODO per cartesian component tau?
            int qfIndex = i * nDimensions + k;
            rNew[i][k] = rOld[i][k] + tau * diffConstant * quantumForceOld[qfIndex] + tau * diffConstant * simpleGaussRandom(idum);
//                    std::cout << "i ndim k " << i * nDimensions + k << " " << quantumForceNew->n_elem << std::endl;
//                    std::cout << "Quantum force " << j << " " << i << " " << k << " " << quantumForceNew[j][i * nDimensions + k] << std::endl;
//                    std::cout << "Old/New: " << rOld[j][i][k] << " " << rNew[j][i][k] << std::endl;
        }
        double ratio = wave->ratio(rNew[i], i);
        // Apply fixed node approximation (keep sign or reject move)
        if(ratio > 0) {
            // Compute weight function
            wave->gradient(rNew, quantumForceNew);
            double argSum = 0;
            for(int j = 0; j < nParticles; j++) {
                for(int k = 0; k < nDimensions; k++) {
                    int qfIndex = j * nDimensions + k;
                    double quantumForceSum = quantumForceOld[qfIndex] + quantumForceNew[qfIndex];
                    double qfPositionDiff = 0.5 * tau * diffConstant * (quantumForceOld[qfIndex] - quantumForceNew[qfIndex]) - (rNew[j][k] - rOld[j][k]);
                    argSum += 0.5 * quantumForceSum * qfPositionDiff;
                }
            }
            double greensRatio = exp(argSum);
            double weight = ratio*ratio * greensRatio;
            // Accept move according to Metropolis probability
            if(weight > ran2(idum)) {
                wave->acceptMove(i);
                rOld[i] = rNew[i];
                quantumForceOld = quantumForceNew;
            } else {
                wave->rejectMove();
                rNew[i] = rOld[i];
            }
        } else {
            wave->rejectMove();
            rNew[i] = rOld[i];
        }
//        wave->initialize(rNew);
        // Compute branching factor PB
        localEnergyNew = hamiltonian->energy(wave, rNew);
//        std::cout << "Energies: " << rNew[0][0] << "," << rNew[0][1] <<" " << rNew[1][0] << "," << rNew[1][1] << ": " << localEnergyNew << " " << localEnergyOld << std::endl;
        double branchingFactor = exp(- tau * (0.5 * (localEnergyOld + localEnergyNew) - trialEnergy));
        // Make int(PB + u) copies
//                std::cout << "PB: " << branchingFactor << " " << trialEnergy << " " << localEnergyOld << " " << localEnergyNew << std::endl;
        int reproductions = (int) (branchingFactor + ran2(idum));
//                std::cout << reproductions << std::endl;
        // Accumulate the energy and any observables weighted by PB
        m_energy += localEnergyNew * branchingFactor;
        m_changeInEnergySamples++;
        if(reproductions == 0) {
//            std::cout << "Killing walker" << std::endl;
            m_aliveNew = false;
        } else {
            if(reproductions > 1) {
                for(int repro = 0; repro < reproductions; repro++) {
//                    std::cout << "Reproduction" << std::endl;
                    bool foundDeadWalker = false;
                    for(int walkerIndex = 0; walkerIndex < nWalkersMax; walkerIndex++) {
                        DiffusionWalker *otherWalker = otherWalkers[walkerIndex];
                        // find a dead walker to ressurect
                        if(!otherWalker->aliveNew() && !foundDeadWalker) {
                            otherWalker->setAliveNew(true);
                            otherWalker->copyOtherWalker(this);
                            foundDeadWalker = true;
                            break;
                        }
                    }
                    if(!foundDeadWalker) {
//                        std::cout << "Could not find an available spot to clone walker" << std::endl;
                    }
                }
            }
        } // END if branching
//                std::cout << "Alive walkers: " << nWalkersAlive << std::endl;
//            std::cout << "Should make " << (int) (branchingFactor + ran2(idum)) << " copies" << std::endl;
        localEnergyOld = localEnergyNew;
    } // END for every particle
}

void DiffusionWalker::progressToNextStep()
{
    Walker::progressToNextStep();
    m_aliveOld = m_aliveNew;
}
