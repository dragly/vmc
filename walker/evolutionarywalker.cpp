#include "evolutionarywalker.h"

#include "../random.h"

EvolutionaryWalker::EvolutionaryWalker(Config *config) :
    Walker(config)
{
}

void EvolutionaryWalker::advance() {
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
            rNew[i][k] = rOld[i][k] + tau * diffConstant * quantumForceOld[qfIndex] + sqrt(2 * tau * diffConstant) * simpleGaussRandom(idum);
//                    std::cout << "i ndim k " << i * nDimensions + k << " " << quantumForceNew->n_elem << std::endl;
//                    std::cout << "Quantum force " << j << " " << i << " " << k << " " << quantumForceNew[j][i * nDimensions + k] << std::endl;
//                    std::cout << "Old/New: " << rOld[j][i][k] << " " << rNew[j][i][k] << std::endl;
        }
        double ratio = wave->ratio(rNew[i], i);
        quantumForceNew *= 2;
        // Apply fixed node approximation (keep sign or reject move)
        if(ratio > 0) {
            // Compute weight function
            wave->gradient(rNew, quantumForceNew);
            quantumForceNew *= 2;
            double argSum = 0;
            for(int k = 0; k < nDimensions; k++) {
                int qfIndex = i * nDimensions + k;
                double quantumForceSum = quantumForceOld[qfIndex] + quantumForceNew[qfIndex];
                double qfPositionDiff = 0.5 * tau * diffConstant * (quantumForceOld[qfIndex] - quantumForceNew[qfIndex]) - (rNew[i][k] - rOld[i][k]);
                argSum += 0.5 * quantumForceSum * qfPositionDiff;
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
        localEnergyOld = hamiltonian->energy(wave, rOld);
        // Accumulate the energy and any observables weighted by PB
        m_energy += localEnergyNew;
        m_changeInEnergySamples++;
//                std::cout << "Alive walkers: " << nWalkersAlive << std::endl;
//            std::cout << "Should make " << (int) (branchingFactor + ran3(idum)) << " copies" << std::endl;
        localEnergyOld = localEnergyNew;
    } // END for every particle
}
