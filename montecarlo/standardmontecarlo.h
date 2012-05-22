#ifndef MONTECARLOSTANDARD_H
#define MONTECARLOSTANDARD_H
#include <armadillo>
using namespace arma;

#include "montecarlo.h"
#include "../hamiltonian/hamiltonian.h"
#include "../wavefunction/wavefunction.h"
#include "../config.h"

class StandardMonteCarlo : public MonteCarlo
{
public:
    StandardMonteCarlo(Config* config);

    void sample(int nCycles);

    ~StandardMonteCarlo();
private:
    int rank;
    bool firstSample;
};

#endif // MONTECARLOSTANDARD_H
