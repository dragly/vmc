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
    void sample(int nCycles, double *energies, double *allEnergies);

    ~MonteCarloMetropolisHastings();
    void quantumForce(vec2 *rPosition, double *forceVectorNew);
private:
    int m_nParticles;
    int m_nDimensions;
    double charge;
    int rank;
    double step_length;
    vec2 *rOld;
    vec2 *rNew;
    double *waveGradientOld;
    double *waveGradientNew;
    double* forceVectorNew;
    double* forceVectorOld;
    double* forceVectorSum;
    double* forceVectorDiff ;
    double* positionDiff ;
    long idum;
    WaveFunction* wave;
    Hamiltonian* hamiltonian;
};

#endif // MONTECARLOMETROPOLISHASTINGS_H
