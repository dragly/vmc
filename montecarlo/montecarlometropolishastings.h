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
private:
    int rank;
    double stepLength;
    vec2 *rOld;
    vec2 *rNew;
//    vec waveGradientOld;
//    vec waveGradientNew;
    vec quantumForceNew;
    vec quantumForceOld;
//    vec forceVectorSum;
//    vec forceVectorDiff ;
    vec2 positionDiff ;
    Hamiltonian* hamiltonian;
    bool firstSample;
};

#endif // MONTECARLOMETROPOLISHASTINGS_H
