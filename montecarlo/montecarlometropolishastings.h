#ifndef MONTECARLOMETROPOLISHASTINGS_H
#define MONTECARLOMETROPOLISHASTINGS_H
#include <armadillo>
using namespace arma;

#include "montecarlo.h"
#include "../hamiltonian/hamiltonian.h"
#include "../wavefunction/wavefunction.h"

class MonteCarloMetropolisHastings : public MonteCarlo
{
public:
    MonteCarloMetropolisHastings(Config *config);
    void sample(int nCycles);

    ~MonteCarloMetropolisHastings();
    void quantumForce(vec2 rPosition[], vec &forceVector);
private:
    int rank;
    double stepLength;
    vec2 *rOld;
    vec2 *rNew;
    vec waveGradientOld;
    vec waveGradientNew;
    vec forceVectorNew;
    vec forceVectorOld;
    vec forceVectorSum;
    vec forceVectorDiff ;
    vec2 positionDiff ;
    WaveFunction* wave;
    Hamiltonian* hamiltonian;
};

#endif // MONTECARLOMETROPOLISHASTINGS_H
