#ifndef MINIMIZEREVOLUTIONARY_H
#define MINIMIZEREVOLUTIONARY_H

#include <vector>
#include <armadillo>
using namespace arma;

#include "minimizer.h"
#include "../evolver/evolver.h"

/*!
    An evolutionary approach to minimizing the energy by selection of
    the parameters. This class is written with the intention to be as
    general as possible with regards to the problem and could hopefully
    be used for a completely different problem.
*/
class MinimizerEvolutionary : public Minimizer, public Evolver
{
public:
    MinimizerEvolutionary(Config *config);

    void runMinimizer();
    void loadConfiguration(INIParser * settings);

    void startEvolution();

private:
    double fitness(vec &coefficients);

    // Monte Carlo stuff
    WaveFunction *wave;
    MonteCarlo *monteCarlo;
    Hamiltonian *hamiltonian;
};

#endif // MINIMIZEREVOLUTIONARY_H
