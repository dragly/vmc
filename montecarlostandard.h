#ifndef MONTECARLOSTANDARD_H
#define MONTECARLOSTANDARD_H

#include "montecarlo.h"
#include "hamiltonian/hamiltonian.h"
#include "wavefunction.h"

class MonteCarloStandard : public MonteCarlo
{
public:
    MonteCarloStandard(WaveFunction* wave, Hamiltonian* hamiltonian, int number_particles, int dimension, double charge, int rank, double step_length);
    void sample(int number_cycles, double *energies);
    double local_energy(double **r, double alpha, double wfold);

    ~MonteCarloStandard();
private:
    int number_particles;
    int dimension;
    double charge;
    int rank;
    double step_length;
    double **r_old, **r_new;
    long idum;
};

#endif // MONTECARLOSTANDARD_H
