#include "evolutionarywalker.h"

EvolutionaryWalker::EvolutionaryWalker(Config *config) :
    Walker(config)
{
}

void EvolutionaryWalker::advance() {
//    for(int i = 0; i < nParticles; i++) {
//        wave->gradient(rOld, i, quantumForceOld);

//        // Propose move (with quantum force)
//        for(int k = 0; k < nDimensions; k++) {
//            // TODO per cartesian component tau?
//            int qfIndex = i * nDimensions + k;
//            rNew[i][k] = rOld[i][k] + tau * diffConstant * quantumForceOld[qfIndex] + tau * diffConstant * simpleGaussRandom(idum);
////                    std::cout << "i ndim k " << i * nDimensions + k << " " << quantumForceNew->n_elem << std::endl;
////                    std::cout << "Quantum force " << j << " " << i << " " << k << " " << quantumForceNew[j][i * nDimensions + k] << std::endl;
////                    std::cout << "Old/New: " << rOld[j][i][k] << " " << rNew[j][i][k] << std::endl;
//        }
//        double ratio = wave->ratio(rNew[i], i);
//        // Apply fixed node approximation (keep sign or reject move)
//        if(ratio > 0) {
//            // Compute weight function

//            wave->gradient(rNew, i, quantumForceNew);
//            double argSum = 0;
//            for(int j = 0; j < nParticles; j++) {
//                for(int k = 0; k < nDimensions; k++) {
//                    int qfIndex = j * nDimensions + k;
//                    double quantumForceSum = quantumForceOld[qfIndex] + quantumForceNew[qfIndex];
//                    double qfPositionDiff = 0.5 * tau * diffConstant * (quantumForceOld[qfIndex] - quantumForceNew[qfIndex]) - (rNew[j][k] - rOld[j][k]);
//                    argSum += 0.5 * quantumForceSum * qfPositionDiff;
//                }
//            }
//            double greensRatio = exp(argSum);
//            double weight = ratio*ratio * greensRatio;
//            // Accept move according to Metropolis probability
//            if(weight > ran2(idum)) {
//                wave->acceptMove(i);
//                rOld[i] = rNew[i];
//            } else {
//                wave->rejectMove();
//                rNew[i] = rOld[i];
//            }
//        } else {
//            wave->rejectMove();
//        }

//        // Compute branching factor PB
//        localEnergyNew = hamiltonian->energy(wave, rNew);
//        //     std::cout << "Alive walkers: " << nWalkersAlive << std::endl;
////            std::cout << "Should make " << (int) (branchingFactor + ran2(idum)) << " copies" << std::endl;
//        localEnergyOld = localEnergyNew;
//    } // END for every particle
}
