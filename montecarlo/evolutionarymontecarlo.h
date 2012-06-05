#ifndef EVOLUTIONARYMONTECARLO_H
#define EVOLUTIONARYMONTECARLO_H

#include "montecarlo.h"
#include "../evolver/evolver.h"
#include "../walker/evolutionarywalker.h"

class Config;

/*!
  * \brief A tested approach to an evolutionary diffusion Monte Carlo algorithm. Doesn't work. Yet.
  */
class EvolutionaryMonteCarlo : public MonteCarlo, public Evolver
{
public:
    EvolutionaryMonteCarlo(Config *config, int nGenes, int nIndividuals, int nPopulations);
    void sample(int nCycles);
    double fitness(vec &genes, int population, int individual);

private:
    int nWalkers;
    EvolutionaryWalker **walkers;
    vec2 *positions;
    int correlationStep;
    double trialEnergy;
    double energySum;
    int nEnergyUpdates;
    vec *meanEnergies;
    vec *meanEnergyCycles;
};

#endif // EVOLUTIONARYMONTECARLO_H
