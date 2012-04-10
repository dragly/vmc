#include "minimizerevolutionary.h"

MinimizerEvolutionary::MinimizerEvolutionary(Config *config) :
    Minimizer(config)
{
}

void MinimizerEvolutionary::runMinimizer()
{
    startEvolution();
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
    // TODO do the monte carlo integral to find the energy for the set coefficients
}

