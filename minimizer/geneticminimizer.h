#ifndef GENETICMINIMIZER_H
#define GENETICMINIMIZER_H

#include <vector>
#include <armadillo>
using namespace arma;

#include "minimizer.h"
#include "../evolver/evolver.h"

class EvolutionaryWalker;

/*!
    An evolutionary approach to minimizing the energy by selection of
    the parameters. This class is written with the intention to be as
    general as possible with regards to the problem and could hopefully
    be used for a completely different problem.
*/
class GeneticMinimizer : public Minimizer, public Evolver
{
public:
    GeneticMinimizer(Config *config);

    void runMinimizer();
    void loadConfiguration(INIParser * settings);

    void setNSamples(int nStart, int nEnd) {
        nSamples = nStart;
        nSamplesStart = nStart;
        nSamplesEnd = nEnd;
    }
    void setNCycles(int nCycles) {
        this->nCycles = nCycles;
    }

private:
    double fitness(vec &coefficients, int population, int individual);


    double nSamples;

    int nSamplesStart;
    int nSamplesEnd;
    int nCycles;
    vec* energies;

    int nParticles;
    int nDimensions;
    double stepLength;
    // Monte Carlo stuff
    MonteCarlo *monteCarlo;
    WaveFunction *wave;
    Hamiltonian *hamiltonian;
};

#endif // GENETICMINIMIZER_H
