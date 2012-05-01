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
    void sample(int nCycles, bool storeEnergies);

    ~MonteCarloMetropolisHastings();
    void quantumForce(vec2 rPosition[], vec2 &forceVector);
private:
    int m_nParticles;
    int m_nDimensions;
    int rank;
    double step_length;
    vec2 *rOld;
    vec2 *rNew;
    vec2 waveGradientOld;
    vec2 waveGradientNew;
    vec2 forceVectorNew;
    vec2 forceVectorOld;
    vec2 forceVectorSum;
    vec2 forceVectorDiff ;
    vec2 positionDiff ;
    WaveFunction* wave;
    Hamiltonian* hamiltonian;
};

#endif // MONTECARLOMETROPOLISHASTINGS_H
