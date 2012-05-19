#include "diffusionmontecarlo.h"
#include "montecarlostandard.h"
#include "../wavefunction/wavefunction.h"
#include "../wavefunction/waveslater.h"
#include "../random.h"
#include "montecarlometropolishastings.h"

DiffusionMonteCarlo::DiffusionMonteCarlo(Config *config) :
    MonteCarlo(config),
    tau(0.01)
{
    nWalkersMax = 10000;
    nWalkersIdeal = 500;
    nWalkersAlive = nWalkersIdeal;
    correlationStep = 200;

    rNew = new vec2*[nWalkersMax];
    for(int i = 0; i < nWalkersMax; i++) {
        rNew[i] = new vec2[nParticles];
    }
    rOld = new vec2*[nWalkersMax];
    for(int i = 0; i < nWalkersMax; i++) {
        rOld[i] = new vec2[nParticles];
    }
    waves = new WaveFunction*[nWalkersMax];
    for(int i = 0; i < nWalkersMax; i++) {
        waves[i] = wave->clone();
    }
    aliveOld = new bool[nWalkersMax];
    aliveNew = new bool[nWalkersMax];
    for(int i = 0; i < nWalkersMax; i++) {
        if(i < nWalkersAlive) {
            aliveOld[i] = true;
            aliveNew[i] = true;
        } else {
            aliveOld[i] = false;
            aliveNew[i] = false;
        }
    }
    std::cout << "Done constructing DMC class" << std::endl;
}

// Note: Implementation from http://www.ornl.gov/~pk7/thesis/pkthnode21.html

void DiffusionMonteCarlo::sample(int nCycles)
{
    // Initialize ensemble of walkers from VMC best guess
    MonteCarloMetropolisHastings *monteCarlo = new MonteCarloMetropolisHastings(config);
    monteCarlo->setRecordMoves(true, nWalkersAlive * nParticles);
    monteCarlo->setTerminalizationEnabled(true);
    monteCarlo->sample(nWalkersAlive * correlationStep);
    double trialEnergy = monteCarlo->energy();
    std::cout << "Initial trial energy was " << trialEnergy << std::endl;

    vec *quantumForceNew = new vec[nWalkersMax];
    vec *quantumForceOld = new vec[nWalkersMax];

    double *localEnergyOld = new double[nWalkersMax];
    double *localEnergyNew = new double[nWalkersMax];

    ofstream scatterfile;
    scatterfile.open("positions-init.dat");
    for(int j = 0; j < nWalkersAlive; j++) {
        for(int i = 0; i < nParticles; i++) {
            vec2 **moves = monteCarlo->moves();
            rOld[j][i] = moves[j][i];
            rNew[j][i] = moves[j][i];
            if(aliveOld[j]) {
                scatterfile << moves[j][i][0] << "\t" << moves[j][i][1] << std::endl;
            }
        }
        quantumForceNew[j] = zeros<vec>(nDimensions * nParticles);
        quantumForceOld[j] = zeros<vec>(nDimensions * nParticles);
        waves[j]->initialize(rNew[j]);
        waves[j]->gradient(rNew[j], 0, quantumForceNew[j]);
        quantumForceOld[j] = quantumForceNew[j];
        localEnergyOld[j] = hamiltonian->energy(waves[j], rNew[j]);
        localEnergyNew[j] = localEnergyOld[j];
    }
    scatterfile.close();

    int blockLength = 100;
    double energySum = 0;
    int nEnergySamples = 0;
    double diffConstant = 0.5;
    // For every cycle:
    for(int cycle = 0; cycle < nCycles; cycle++) {
        // For every walker (configuration)
        for(int walker = 0; walker < nWalkersMax; walker++) {
            if(aliveOld[walker]) {
                // For every particle
                for(int i = 0; i < nParticles; i++) {
                    waves[walker]->gradient(rOld[walker], i, quantumForceOld[walker]);

                    // Propose move (with quantum force)
                    for(int k = 0; k < nDimensions; k++) {
                        // TODO per cartesian component tau?
                        int qfIndex = i * nDimensions + k;
                        rNew[walker][i][k] = rOld[walker][i][k] + tau * diffConstant * quantumForceOld[walker][qfIndex] + tau * diffConstant * simpleGaussRandom(idum);
    //                    std::cout << "i ndim k " << i * nDimensions + k << " " << quantumForceNew->n_elem << std::endl;
    //                    std::cout << "Quantum force " << j << " " << i << " " << k << " " << quantumForceNew[j][i * nDimensions + k] << std::endl;
    //                    std::cout << "Old/New: " << rOld[j][i][k] << " " << rNew[j][i][k] << std::endl;
                    }
                    double ratio = waves[walker]->ratio(rNew[walker][i], i);
                    // Apply fixed node approximation (keep sign or reject move)
                    if(ratio > 0) {
                        // Compute weight function

                        waves[walker]->gradient(rNew[walker], i, quantumForceNew[walker]);
                        double argSum = 0;
                        for(int j = 0; j < nParticles; j++) {
                            for(int k = 0; k < nDimensions; k++) {
                                int qfIndex = j * nDimensions + k;
                                double quantumForceSum = quantumForceOld[walker][qfIndex] + quantumForceNew[walker][qfIndex];
                                double qfPositionDiff = 0.5 * tau * diffConstant * (quantumForceOld[walker][qfIndex] - quantumForceNew[walker][qfIndex]) - (rNew[walker][j][k] - rOld[walker][j][k]);
                                argSum += 0.5 * quantumForceSum * qfPositionDiff;
                            }
                        }
                        double greensRatio = exp(argSum);
                        double weight = ratio*ratio * greensRatio;
                        // Accept move according to Metropolis probability
                        if(weight > ran2(idum)) {
                            waves[walker]->acceptMove(i);
                            rOld[walker][i] = rNew[walker][i];
                        } else {
                            waves[walker]->rejectMove();
                            rNew[walker][i] = rOld[walker][i];
                        }
                    } else {
                        waves[walker]->rejectMove();
                    }

                    // Compute branching factor PB
                    localEnergyNew[walker] = hamiltonian->energy(waves[walker], rNew[walker]);
                    double branchingFactor = exp(- tau * (0.5 * (localEnergyOld[walker] + localEnergyNew[walker]) - trialEnergy));
                    // Make int(PB + u) copies
    //                std::cout << "PB: " << branchingFactor << " " << trialEnergy << " " << localEnergyOld[walker] << " " << localEnergyNew[walker] << std::endl;
                    int reproductions = (int) (branchingFactor + ran2(idum));
    //                std::cout << reproductions << std::endl;
                    if(reproductions == 0) {
    //                    std::cout << "Killing walker" << std::endl;
                        nWalkersAlive--;
                        aliveNew[walker] = false;
                    } else {
                        // Accumulate the energy and any observables weighted by PB
                        energySum += localEnergyNew[walker] * branchingFactor;
                        nEnergySamples++;
                        if(reproductions > 1) {
                            for(int repro = 0; repro < reproductions; repro++) {
    //                            std::cout << "Good choice!" << std::endl;
                                bool foundDeadWalker = false;
                                for(int walker2 = 0; walker2 < nWalkersMax; walker2++) {
                                    // find a dead walker to ressurect
                                    if(!aliveNew[walker2] && !foundDeadWalker) {
                                        nWalkersAlive++;
                                        aliveNew[walker2] = true;
                                        foundDeadWalker = true;
                                        for(int i = 0; i < nParticles; i++) {
                                            rNew[walker2][i] = rNew[walker][i];
                                            rOld[walker2][i] = rOld[walker][i];
                                        }
                                        quantumForceNew[walker2] = quantumForceNew[walker];
                                        quantumForceOld[walker2] = quantumForceOld[walker];
                                        waves[walker2]->initialize(rNew[walker2]);
                                        localEnergyNew[walker2] = localEnergyNew[walker];
                                        localEnergyOld[walker2] = localEnergyOld[walker];

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
                    localEnergyOld[walker] = localEnergyNew[walker];
                } // END for every particle
            } // END if walker alive
        } // END for each walker
//        std::cout << "Alive walkers: " << nWalkersAlive << "\xd" << std::endl;
        // Repeat configuration moves for about 100 - 1000 steps
        if(cycle < 400 || !(cycle % blockLength)) {
            // Update trial energy ET to bring it closer to the current ensemble
            trialEnergy = energySum / nEnergySamples;
            std::cout << "Trial energy is now " << trialEnergy << " with " << nWalkersAlive << " walkers at cycle " << cycle << std::endl;
            // Renormalise the number of walkers to the target number by creating or deleting walkers
//            while(nWalkersAlive > nWalkersIdeal) {
//                int randomWalker = ran2(idum) * nWalkersMax;
//                if(aliveNew[randomWalker]) {
//                    aliveNew[randomWalker] = false;
//                    nWalkersAlive--;
//                }
//            }

            energySum = 0;
            nEnergySamples = 0;
        }
        for(int walker = 0; walker < nWalkersMax; walker++) {
//            std::cout << "Walkers alive " << aliveOld[walker] << " " << aliveNew[walker] << std::endl;
            aliveOld[walker] = aliveNew[walker];
//            std::cout << "Walkers alive after " << aliveOld[walker] << " " << aliveNew[walker] << std::endl;
        }
        // Repeat and rinse
    }
    scatterfile.open("positions-end.dat");
    for(int j = 0; j < nWalkersMax; j++) {
        for(int i = 0; i < nParticles; i++) {
            if(aliveOld[j]) {
                scatterfile << rNew[j][i][0] << "\t" << rNew[j][i][1] << std::endl;
            }
        }
    }
    scatterfile.close();

    m_energy = trialEnergy;
}
