#ifndef MONTECARLOSTANDARD_H
#define MONTECARLOSTANDARD_H
#include <armadillo>
using namespace arma;

#include "montecarlo.h"
#include "../hamiltonian/hamiltonian.h"
#include "../wavefunction/wavefunction.h"
#include "../config.h"

class MonteCarloStandard : public MonteCarlo
{
public:
    MonteCarloStandard(Config* config);

    void sample(int nCycles);

    ~MonteCarloStandard();
private:
    int rank;
    double stepLength;
    vec2 *rOld;
    vec2 *rNew;
    bool firstSample;
};

#endif // MONTECARLOSTANDARD_H
