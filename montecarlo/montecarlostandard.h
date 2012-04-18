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
    void sample(int nCycles, double *energies, double *allEnergies);

    ~MonteCarloStandard();
private:
    int m_nParticles;
    int m_nDimensions;
    double charge;
    int rank;
    double step_length;
    vec2 *rOld;
    vec2 *rNew;
    long idum;
};

#endif // MONTECARLOSTANDARD_H
