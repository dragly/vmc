#include "evolutionarymontecarlo.h"
#include "../montecarlo/montecarlometropolishastings.h"
#include "../config.h"
#include "../random.h"

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
}

double EvolutionaryMonteCarlo::fitness(vec &genes, int population, int individual) {
    int nWalkerSamples = 5000;
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
    for(int sample = 0; sample < nWalkerSamples; sample++) {
        for(int i = 0; i < nWalkers; i++) {
            EvolutionaryWalker* walker = walkers[i];
            walker->advance();
            totalEnergy += walker->energy();
            nTotalSamples += walker->changeInEnergySamples();
        }
    }
    double meanEnergy = totalEnergy / (nTotalSamples);
    meanEnergies[population][individual] = meanEnergy;
    // return difference between trial energy and this individuals energy
    return fabs(meanEnergy - trialEnergy);
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
            }
        }
    }
    for(int cycle = 0; cycle < nCycles; cycle++) {
        energySum = 0;
        nEnergyUpdates = 0;
        allBestValue = INFINITY;
        evolve(10, 100);
        std::cout << "All best value was " << allBestValue << " found in " << meanEnergies[allBestPopulationIndex][allBestIndex] << std::endl;
        for(int i = 0; i < nPopulations; i++) {
            // Avoid the random newcomers by selecting the first three quarters (the best and their children)
            for(uint j = 0; j < bestIndices[i].size() * 3. / 4.; j++) {
                int index = bestIndices[i][j];
                std::cout << "Mean energies: " << meanEnergies[i][index] << std::endl;
                energySum += meanEnergies[i][index];
                nEnergyUpdates++;
            }
        }
        trialEnergy = energySum / nEnergyUpdates;
        std::cout << "New trial energy " << trialEnergy << std::endl;
    }
    m_energy = trialEnergy;
}
