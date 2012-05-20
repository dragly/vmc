#include "evolutionarymontecarlo.h"
#include "../config.h"

EvolutionaryMonteCarlo::EvolutionaryMonteCarlo(Config *config, int nGenes, int nIndividuals, int nPopulations) :
    MonteCarlo(config),
    Evolver(nGenes, nIndividuals, nPopulations)
{

    nWalkers = nGenes / nParticles;
}

double EvolutionaryMonteCarlo::fitness(vec &coefficients) {
    // take all the coefficients, spawn particles in these positions
    for(int i = 0; i < nWalkers; i++) {
        // walker
    }
    // sample energies around these positions
    // return difference between trial energy and this individuals energy
    return 0;
}

void EvolutionaryMonteCarlo::sample(int nCycles)
{
    evolve(nCycles, 100);
}
