#ifndef MONTECARLOMETROPOLISHASTINGS_H
#define MONTECARLOMETROPOLISHASTINGS_H

#include "montecarlo.h"
#include "../hamiltonian/hamiltonian.h"
#include "../wavefunction.h"

class MonteCarloMetropolisHastings : public MonteCarlo
{
public:
    MonteCarloMetropolisHastings(WaveFunction* m_wave, Hamiltonian* m_hamiltonian, int number_particles, int dimension, double charge, int rank, double step_length);
    void sample(int nCycles, double *energies, double *allEnergies);

    ~MonteCarloMetropolisHastings();
private:
    int number_particles;
    int dimension;
    double charge;
    int rank;
    double step_length;
    double **r_old;
    double **r_new;
    long idum;
};

#endif // MONTECARLOMETROPOLISHASTINGS_H
