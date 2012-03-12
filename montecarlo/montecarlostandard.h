#ifndef MONTECARLOSTANDARD_H
#define MONTECARLOSTANDARD_H

#include "montecarlo.h"
#include "../hamiltonian/hamiltonian.h"
#include "../wavefunction.h"

class MonteCarloStandard : public MonteCarlo
{
public:
    MonteCarloStandard(WaveFunction* m_wave, Hamiltonian* m_hamiltonian, int number_particles, int dimension, double charge, int rank, double step_length);
    void sample(int nCycles, double *energies, double *allEnergies);

    ~MonteCarloStandard();
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

#endif // MONTECARLOSTANDARD_H
