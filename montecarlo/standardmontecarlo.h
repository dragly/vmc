#ifndef MONTECARLOSTANDARD_H
#define MONTECARLOSTANDARD_H
#include <armadillo>
using namespace arma;

#include "montecarlo.h"
#include "../hamiltonian/hamiltonian.h"
#include "../wavefunction/wavefunction.h"
#include "../config.h"

/*!
  * \brief Implements the brute force Monte Carlo algorithm.
  */
class StandardMonteCarlo : public MonteCarlo
{
public:
    StandardMonteCarlo(Config* config);

    void sample(int nSamples);

    ~StandardMonteCarlo();
private:
    int myRank;
    bool firstSample;
};

#endif // MONTECARLOSTANDARD_H
