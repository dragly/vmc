#include "diffusionmontecarlo.h"
#include "montecarlostandard.h"
#include "../wavefunction/wavefunction.h"
#include "../wavefunction/waveslater.h"
#include "../random.h"

DiffusionMonteCarlo::DiffusionMonteCarlo(Config *config) :
    MonteCarlo(config),
    tau(1)
{
    nWalkersMax = 10000;
    nWalkersIdeal = 100;
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
    MonteCarloStandard *monteCarloStandard = new MonteCarloStandard(config);
    monteCarloStandard->setRecordMoves(true, nWalkersMax * nParticles);
    monteCarloStandard->sample(nWalkersMax * correlationStep);
    double trialEnergy = monteCarloStandard->energy();

    vec *quantumForceNew = new vec[nWalkersMax];
    vec *quantumForceOld = new vec[nWalkersMax];

    double *localEnergyOld = new double[nWalkersMax];
    double *localEnergyNew = new double[nWalkersMax];

    ofstream scatterfile;
    scatterfile.open("positions-init.dat");
    for(int j = 0; j < nWalkersMax; j++) {
        for(int i = 0; i < nParticles; i++) {
            vec2 **moves = monteCarloStandard->moves();
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

    // For every cycle:
    for(int cycle = 0; cycle < nCycles; cycle++) {
        // For every walker (configuration)
        for(int walker = 0; walker < nWalkersMax; walker++) {
            if(aliveOld[walker]) {
                // For every particle
                for(int i = 0; i < nParticles; i++) {
                    waves[walker]->gradient(rNew[walker], i, quantumForceNew[walker]);

                    // Propose move (with quantum force)
                    for(int k = 0; k < nDimensions; k++) {
                        // TODO per cartesian component tau?
                        rNew[walker][i][k] = rOld[walker][i][k] + tau * quantumForceNew[walker][i * nDimensions + k] + simpleGaussRandom(idum);
    //                    std::cout << "i ndim k " << i * nDimensions + k << " " << quantumForceNew->n_elem << std::endl;
    //                    std::cout << "Quantum force " << j << " " << i << " " << k << " " << quantumForceNew[j][i * nDimensions + k] << std::endl;
    //                    std::cout << "Old/New: " << rOld[j][i][k] << " " << rNew[j][i][k] << std::endl;
                    }
                    double ratio = waves[walker]->ratio(rNew[walker][i], i);
                    // Apply fixed node approximation (keep sign or reject move)
                    if(ratio > 0) {
                        // Compute weight function
                        double argSum = 0;
                        for(int k = 0; k < nDimensions; k++) {
                            double argument = rNew[walker][i][k] - rOld[walker][i][k] - tau * quantumForceOld[walker][i * nDimensions + k];
                            argSum += argument * argument;
                        }
                        argSum *= 1 / 4.;
                        double greensRatio = pow((2 * M_PI * tau),(-3*nParticles / 2)) * exp(-(argSum) / 2 * tau);
                        double weight = ratio*ratio * greensRatio;
                        // Accept move according to Metropolis probability
                        if(weight > ran2(idum) - 0.5) {
                            waves[walker]->acceptMove(i);
                            rOld[walker][i] = rNew[walker][i];
                        } else {
                            waves[walker]->rejectMove();
                            rNew[walker][i] = rOld[walker][i];
                        }
                    } else {
                        waves[walker]->rejectMove();
                    }
                } // END for every particle
                // Compute branching factor PB
                localEnergyNew[walker] = hamiltonian->energy(waves[walker], rNew[walker]);
                double branchingFactor = exp(- tau * (1 / 2. * (localEnergyOld[walker] + localEnergyNew[walker]) - trialEnergy));
                // Accumulate the energy and any observables weighted by PB
                // Make int(PB + u) copies
                if(branchingFactor + ran2(idum) > 1) {
//                    std::cout << "Good choice!" << std::endl;
                    bool foundDeadWalker = false;
                    for(int walker2 = 0; walker2 < nWalkersMax; walker2++) {
                        // find a dead walker to ressurect
                        if(!aliveNew[walker2]) {
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
                        std::cout << "Could not find an available spot to clone walker" << std::endl;
                    }
                } else {
//                    std::cout << "Didn't reproduce" << std::endl;
                    nWalkersAlive--;
                    aliveNew[walker] = false;
                } // END if branching
//                std::cout << "Alive walkers: " << nWalkersAlive << std::endl;
    //            std::cout << "Should make " << (int) (branchingFactor + ran2(idum)) << " copies" << std::endl;
                localEnergyOld[walker] = localEnergyNew[walker];
            } // END if walker alive
        } // END for each walker
        // Repeat configuration moves for about 100 - 1000 steps
        if(!(cycle % 10)) {
            // Update trial energy ET to bring it closer to the current ensemble
            double energySum = 0;
            for(int walker = 0; walker < nWalkersMax; walker++) {
                energySum += localEnergyNew[walker];
            }
            trialEnergy = energySum / nWalkersMax;
            std::cout << "Trial energy is now " << trialEnergy << std::endl;
            std::cout << "Alive walkers: " << nWalkersAlive << std::endl;
            // Renormalise the number of walkers to the target number by creating or deleting walkers
            while(nWalkersAlive > nWalkersIdeal) {
                int randomWalker = ran2(idum) * nWalkersMax;
                if(aliveNew[randomWalker]) {
                    aliveNew[randomWalker] = false;
                    nWalkersAlive--;
                }
            }
        }
        for(int walker = 0; walker < nWalkersMax; walker++) {
            aliveOld[walker] = aliveNew[walker];
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

    double energySum = 0;
    for(int walker = 0; walker < nWalkersMax; walker++) {
        energySum += localEnergyNew[walker];
    }
    m_energy = energySum / nWalkersMax;
}
