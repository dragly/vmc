#include "minimizerevolutionary.h"
#include "../inih/ini.h"

MinimizerEvolutionary::MinimizerEvolutionary(Config *config) :
    Minimizer(config)
{
    monteCarlo = MonteCarlo::fromName(config->monteCarloClass(), config);
    wave = config->wave();
    hamiltonian = config->hamiltonian();
}

void MinimizerEvolutionary::runMinimizer()
{
    setScaleLimits(-1, 1);
    evolve(1000,250);
}

void MinimizerEvolutionary::loadConfiguration(INIParser *settings)
{
    nIndividuals = settings->GetInteger("MinimizerEvolutionary","nIndividuals",64);
    nPopulations = settings->GetInteger("MinimizerEvolutionary","nPopulations",2);
    nGenes = settings->GetInteger("MinimizerEvolutionary","nGenes",2);
    constructor(nGenes, nPopulations, nIndividuals);
    wave = config->wave();
    hamiltonian = config->hamiltonian();
}

void MinimizerEvolutionary::startEvolution()
{
}

double MinimizerEvolutionary::fitness(vec &coefficients)
{
    double parameters[2];
    parameters[0] = coefficients[0];
    parameters[1] = coefficients[1];
//    std::cout << "Testing parameters " << parameters[0] << " " << parameters[1] << std::endl;
    wave->setParameters(parameters);
    monteCarlo->setTerminalizationEnabled(true);
    monteCarlo->sample(1000);
    // update the energy average and its squared
    double energy = monteCarlo->energy();
//    std::cout << energy << std::endl;
//    cumulativeEnergySquared(i,j) = m_monteCarlo->energySquared();
//    std::cout << "Got energy of " << cumulativeEnergy(i,j) << std::endl;
//    parameter0Map(i,j) = parameters[0];
//    parameter1Map(i,j) = parameters[1];
//    parameters[1] += betaStep;
    return energy;
}

