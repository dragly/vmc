#ifndef MINIMIZEREVOLUTIONARY_H
#define MINIMIZEREVOLUTIONARY_H

#include <vector>
#include <armadillo>
using namespace arma;

#include "minimizer.h"

/*!
    An evolutionary approach to minimizing the energy by selection of
    the parameters. This class is written with the intention to be as
    general as possible with regards to the problem and could hopefully
    be used for a completely different problem.
*/
class MinimizerEvolutionary : Minimizer
{
public:
    MinimizerEvolutionary(Config *config);

    void runMinimizer();
    void loadConfiguration(INIParser * settings);

    void startEvolution();

private:
    double value(vec *coefficients);

    // defines the number of individuals in a population
    // and the size of the genome of each individual
    // read: The number of test functions and coefficients
    double nIndividuals;
    double nGenes;
    double nPopulations;

    // Keep track of the number of cycles
    int currentCycle;
    int nCycles;


    // defines the scaling of the coefficients
    // and how this is limited. Could be coefficients from 1 to 1e6
    // or from 1e3 to 1e4 for example
    double scale;
    double lowScaleLimit;
    double highScaleLimit;
    // at what precision should we start calling for rescales?
    double rescalePrecisionLimit;

    // holds the populations, which each is a vector of genes
    // read: the functions, holding their coefficients
    vec **populations;

    // holds the scoreboard information
    // which individual was best within each population?
    int *bestIndex;
    double *bestValue;
};

#endif // MINIMIZEREVOLUTIONARY_H
