#ifndef EVOLVER_H
#define EVOLVER_H

#include <armadillo>

using namespace arma;

class Evolver
{
public:
    Evolver();
    Evolver(int nGenes, int nIndividuals, int nPopulations);

    void constructor(int nGenes, int nIndividuals, int nPopulations);

    void evolve(int nSteps, int populationMatching);

    void setScaleLimits(double low, double high) {
        lowScaleLimit = low;
        highScaleLimit = high;
    }

    vec allBestGenes;

    ~Evolver();
    void rescale();
protected:
    virtual double fitness(vec &coefficients) = 0;
    void updateBest();

    // defines the number of individuals in a population
    // and the size of the genome of each individual
    // read: The number of test functions and coefficients
    int nIndividuals;
    int nGenes;
    int nPopulations;

    // Keep track of the number of cycles
    int currentCycle;

    long idum;

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
    // which individuals are the best within each population
    uvec *bestIndices;
    vec *values;

    double allBestValue;
    int allBestIndex;
    int allBestPopulationIndex;

    // holds the number of cycles since last improvement
    int cyclesSinceLastImprovement;
    int cyclesSinceLastRescale;

    double lastWorkingScale;
};

#endif // EVOLVER_H
