#include "minimizerevolutionary.h"
#include "../inih/ini.h"

MinimizerEvolutionary::MinimizerEvolutionary(Config *config) :
    Minimizer(config)
{

}

void MinimizerEvolutionary::runMinimizer()
{
    startEvolution();
}

void MinimizerEvolutionary::loadConfiguration(ini *settings)
{
    nIndividuals = settings->GetInteger("MinimizerEvolutionary","nIndividuals",64);
    nPopulations = settings->GetInteger("MinimizerEvolutionary","nPopulations",2);
    nGenes = settings->GetInteger("MinimizerEvolutionary","nGenes",32);
}

void MinimizerEvolutionary::startEvolution()
{
}

/*!
  Calculates everything and returns the current value of the selected
  configuration. Read: The energy for the selected coefficients.
  */
double MinimizerEvolutionary::value(vec *coefficients)
{
    (void)coefficients;
    // TODO do the monte carlo integral to find the energy for the set coefficients
    return 0.0;
}

