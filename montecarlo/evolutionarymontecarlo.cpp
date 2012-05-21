#include "evolutionarymontecarlo.h"
#include "../montecarlo/montecarlometropolishastings.h"
#include "../config.h"
#include "../random.h"

#include <assert.h>

EvolutionaryMonteCarlo::EvolutionaryMonteCarlo(Config *config, int nGenes_, int nIndividuals, int nPopulations) :
    MonteCarlo(config),
    Evolver(nGenes_, nIndividuals, nPopulations)
{
    if(nGenes_ % (nParticles * nDimensions)) {
        std::cerr << "The number of genes must be a multiplum of (nParticles * nDimensions)" << std::endl;
        exit(923);
    }
    correlationStep = 200;
    nWalkers = nGenes_ / (nParticles * nDimensions);
    walkers = new EvolutionaryWalker*[nWalkers];
    for(int i = 0; i < nWalkers; i++) {
        walkers[i] = new EvolutionaryWalker(config);
    }
    positions = new vec2[nParticles];
    meanEnergies = new vec[nPopulations];
    for(int i = 0; i < nPopulations; i++) {
        meanEnergies[i] = zeros<vec>(nIndividuals);
    }
    meanEnergyCycles = new vec[nPopulations];
    for(int i = 0; i < nPopulations; i++) {
        meanEnergyCycles[i] = zeros<vec>(nIndividuals);
    }
}

double EvolutionaryMonteCarlo::fitness(vec &genes, int population, int individual) {
    int nWalkerSamples = 200;
    // take all the coefficients, spawn particles in these positions
    for(int i = 0; i < nWalkers; i++) {
        for(int j = 0; j < nParticles; j++) {
            for(int k = 0; k < nDimensions; k++) {
                int geneIndex = i * nParticles * nDimensions + j * nDimensions + k;
                positions[j][k] = genes[geneIndex];
            }
        }
        walkers[i]->initialize(positions);
    }
    // sample energies around these positions
    double totalEnergy = 0;
    int nTotalSamples = 0;
    double branchingSum = 0;
    for(int sample = 0; sample < nWalkerSamples; sample++) {
        for(int i = 0; i < nWalkers; i++) {

            EvolutionaryWalker* walker = walkers[i];
            double oldEnergy = walker->energy();
            walker->advance();
//            if(sample > nWalkerSamples / 10) {
                double newEnergy = walker->energy();
                if(!isnan(newEnergy)) {
                    totalEnergy += walker->energy();
                    nTotalSamples += walker->changeInEnergySamples();
                }

                double branchingFactor = -exp(- 0.001 * (0.5 * (oldEnergy + newEnergy) - trialEnergy));
                branchingSum += (int)( branchingFactor + ran2(idum));
//            }
        }
    }
    for(int i = 0; i < nWalkers; i++) {
        walkers[i]->advance();
    }
    // Set the new genes to the new positions (think of this as an individual keeping it's new properties after some time)
    for(int k = 0; k < nGenes; k++) {
        int walkerIndex = k / (nParticles * nDimensions);
        int kMarked = (k - walkerIndex * (nDimensions * nParticles)); // kMarked is now between 0 and nDimensions * nParticles
        int particleIndex = kMarked / nDimensions;
        int positionIndex = k % nDimensions;
        vec2 *positionsNew = walkers[walkerIndex]->positionsNew();
        populations[population][individual][k] = positionsNew[particleIndex][positionIndex];
    }
    double meanEnergy = totalEnergy / (nTotalSamples);
    if(nTotalSamples > 0) {
//        meanEnergies[population][individual] += meanEnergy;
        meanEnergies[population][individual] = meanEnergy;
        meanEnergyCycles[population][individual] += 1;
    } else {
        meanEnergy = INFINITY;
    }

//    return exp(- 0.1 * (fabs(meanEnergy - trialEnergy)));
    // return difference between trial energy and this individuals energy
    return fabs(meanEnergy - trialEnergy) * (1 - 0.5 *ran2(idum));
//    return -branchingSum;
}

void EvolutionaryMonteCarlo::sample(int nCycles)
{
    // Initialize ensemble of walkers from VMC best guess
    MonteCarloMetropolisHastings *initialMonteCarlo = new MonteCarloMetropolisHastings(config);
    initialMonteCarlo->setRecordMoves(true, nPopulations * nIndividuals * nParticles * nWalkers);
    initialMonteCarlo->setTerminalizationEnabled(true);
    initialMonteCarlo->sample(nWalkers * nIndividuals * nPopulations * correlationStep);
    trialEnergy = initialMonteCarlo->energy();
    std::cout << "Initial trial energy was " << trialEnergy << std::endl;

    vec2 **moves = initialMonteCarlo->moves();

    // Initialize genes with good guesses
    for(int i = 0; i < nPopulations; i++) {
        for(int j = 0; j < nIndividuals; j++) {
            for(int k = 0; k < nGenes; k++) {
                int walkerIndex = k / (nParticles * nDimensions);
                int kMarked = (k - walkerIndex * (nDimensions * nParticles)); // kMarked is now between 0 and nDimensions * nParticles
                int particleIndex = kMarked / nDimensions;
                int positionIndex = k % nDimensions;
                int moveIndex = walkerIndex + nWalkers * j;
                populations[i][j][k] = moves[moveIndex][particleIndex][positionIndex];
//                populations[i][j][k] = ran2(idum) * 8;
            }
        }
    }
//    trialEnergy = 2;
    for(int cycle = 0; cycle < nCycles; cycle++) {
        allBestValue = INFINITY;
        evolve(10, 100);
        std::cout << "All best value was " << allBestValue << " found in " << meanEnergies[allBestPopulationIndex][allBestIndex] << std::endl;
//        if(allBestValue < 1e-3 || cycle < 5) {
            energySum = 0;
            nEnergyUpdates = 0;
            for(int i = 0; i < nPopulations; i++) {
                // Avoid the random newcomers by selecting the first three quarters (the best and their children)
                for(uint j = 0; j < bestIndices[i].size() * 2. / 4.; j++) {
                    int index = bestIndices[i][j];
//                    int index = j;
                    double meanEnergy = meanEnergies[i][index] /*/ meanEnergyCycles[i][index]*/;
//                    std::cout << "Mean energies: " << index << "\t" << meanEnergy <<  "\tdiff: " << values[i][index] << std::endl;
                    energySum += meanEnergy;
                    nEnergyUpdates++;
                    meanEnergies[i][index] = 0;
                    meanEnergyCycles[i][index] = 0;
                }
            }
            trialEnergy = energySum / nEnergyUpdates;
//            trialEnergy = lowestEnergy;
            std::cout << "New trial energy " << trialEnergy << std::endl;
//        }
    }
    m_energy = trialEnergy;
}
