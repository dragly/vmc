#ifndef MONTECARLOMETROPOLISHASTINGS_H
#define MONTECARLOMETROPOLISHASTINGS_H

#include "montecarlo.h"
#include "../hamiltonian/hamiltonian.h"
#include "../wavefunction.h"

class MonteCarloMetropolisHastings : public MonteCarlo
{
public:
    MonteCarloMetropolisHastings(WaveFunction* wave, Hamiltonian* hamiltonian, int nParticles, int nDimensions, double charge, int rank, double step_length);
    void sample(int nCycles, double *energies, double *allEnergies);

    ~MonteCarloMetropolisHastings();
    void quantumForce(double **rPosition, double *forceVectorNew);
private:
    int m_nParticles;
    int m_nDimensions;
    double charge;
    int rank;
    double step_length;
    double **rOld;
    double **rNew;
    double *waveGradientOld;
    double *waveGradientNew;
    double* forceVectorNew;
    double* forceVectorOld;
    double* forceVectorSum;
    double* forceVectorDiff ;
    double* positionDiff ;
    long idum;
};

#endif // MONTECARLOMETROPOLISHASTINGS_H
