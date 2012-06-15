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

/*!
 * \brief DiffusionWalker::advance suggests a move, accepts/rejects it, samples the energy, kills itself or spawns walkers if necessary
 * \param trialEnergy
 */
void DiffusionWalker::advance(double trialEnergy) {
    Walker::advance();
    m_changeInEnergySamples = 0;
    m_energy = 0;
    m_acceptances = 0;
    m_rejections = 0;
    for(int i = 0; i < nParticles && m_aliveNew; i++) {
//        wave->gradient(rOld, quantumForceOld);
//        quantumForceOld *= 2;
        // Propose move (with quantum force)
        for(int k = 0; k < nDimensions; k++) {
            // TODO per cartesian component tau?
            int qfIndex = i * nDimensions + k;
//            rNew[i][k] = rOld[i][k] + timeStep * diffConstant * quantumForceOld[qfIndex] + sqrt(2 * diffConstant * timeStep) * simpleGaussRandom(idum);
            rNew[i][k] = rOld[i][k] + timeStep * diffConstant * quantumForceOld[qfIndex] + 2 * diffConstant * sqrt(timeStep) * gaussianDeviate(idum);

        }
        wave->prepareGradient(rNew[i], i);
        wave->gradient(rNew, quantumForceNew);
        quantumForceNew *= 2;
        double ratio = wave->ratio(rNew[i], i);
        // Apply fixed node approximation (keep sign or reject move)
        if(ratio > 0) {
            // Compute weight function
            double argSum = 0;
            for(int j = 0; j < nParticles; j++) { // TODO figure out if it is necessary to do this with all particles (is it not zero for non-moved particles?).
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
//                    double squareArgNew = (rOld[j][k] - rNew[j][k] - timeStep * quantumForceNew[qfIndex] / 2);
//                    double newGreen = exp(-squareArgNew*squareArgNew / (2*timeStep));
//                    double squareArgOld = (rNew[j][k] - rOld[j][k] - timeStep * quantumForceOld[qfIndex] / 2);
//                    double oldGreen = exp(-squareArgOld*squareArgOld / (2*timeStep));
//                    greensRatio *= newGreen / oldGreen;
//                }
//            }

            double weight = ratio*ratio * greensRatio;
            // Accept move according to Metropolis probability
            if(ran3(idum) < weight) {
                rOld[i] = rNew[i];
                wave->acceptMove(i);
                quantumForceOld = quantumForceNew;
                m_acceptances++;
                localEnergyNew = hamiltonian->energy(wave, rOld);
            } else {
                wave->rejectMove();
                rNew[i] = rOld[i];
                m_rejections++;
                quantumForceNew = quantumForceOld;
                localEnergyNew = localEnergyOld;
            }
        } else {
            wave->rejectMove();
            rNew[i] = rOld[i];
            quantumForceNew = quantumForceOld;
            localEnergyNew = localEnergyOld;
        }

        // TODO: Try to remove the branching and have a look at the distribution afterwards
        // Compute branching factor PB
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
                for(int repro = 1; repro < reproductions; repro++) {
                    for(int walkerIndex = 0; walkerIndex < nWalkersMax; walkerIndex++) {
                        DiffusionWalker *ressurectedWalker = otherWalkers[walkerIndex];
                        // find a dead walker to ressurect
                        if(!ressurectedWalker->aliveNew()) {
                            ressurectedWalker->setAliveNew(true);
                            ressurectedWalker->copyOtherWalker(this);
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
