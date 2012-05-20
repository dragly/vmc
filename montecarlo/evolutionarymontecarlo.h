#ifndef EVOLUTIONARYMONTECARLO_H
#define EVOLUTIONARYMONTECARLO_H

#include "montecarlo.h"
#include "../evolver/evolver.h"

class Config;

class EvolutionaryMonteCarlo : public MonteCarlo, public Evolver
{
public:
    EvolutionaryMonteCarlo(Config *config, int nGenes, int nIndividuals, int nPopulations);
    void sample(int nCycles);
    double fitness(vec &coefficients);

private:
    int nWalkers;
};

#endif // EVOLUTIONARYMONTECARLO_H
