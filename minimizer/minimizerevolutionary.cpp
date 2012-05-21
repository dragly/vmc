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
    for(int cycle = 0; cycle < 1000; cycle++) {
        allBestValue = INFINITY;
        evolve(1,5);
    }
}

void MinimizerEvolutionary::loadConfiguration(INIParser *settings)
{
    nIndividuals = settings->GetInteger("MinimizerEvolutionary","nIndividuals",16);
    nPopulations = settings->GetInteger("MinimizerEvolutionary","nPopulations",2);
    nGenes = settings->GetInteger("MinimizerEvolutionary","nGenes",2);
    constructor(nGenes, nPopulations, nIndividuals);
    wave = config->wave();
    hamiltonian = config->hamiltonian();
}

void MinimizerEvolutionary::startEvolution()
{
}

double MinimizerEvolutionary::fitness(vec &coefficients, int population, int individual)
{
    double parameters[2];
    parameters[0] = coefficients[0];
    parameters[1] = coefficients[1];
    wave->setParameters(parameters);
    monteCarlo->setTerminalizationEnabled(true);
    monteCarlo->sample(800000);
    // update the energy average and its squared
    double energy = monteCarlo->energy();
    std::cout << "Population\t" << population << ", individual\t" << individual << ", testing parameters\t" << parameters[0] << "\t" << parameters[1] << ", got energy\t" << energy << std::endl;
//    std::cout << energy << std::endl;
//    cumulativeEnergySquared(i,j) = m_monteCarlo->energySquared();
//    std::cout << "Got energy of " << cumulativeEnergy(i,j) << std::endl;
//    parameter0Map(i,j) = parameters[0];
//    parameter1Map(i,j) = parameters[1];
//    parameters[1] += betaStep;
    if(energy == 0) {
        return INFINITY;
    }
    return energy;
}

