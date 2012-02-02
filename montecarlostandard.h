#ifndef MONTECARLOSTANDARD_H
#define MONTECARLOSTANDARD_H

#include "montecarlo.h"

class MonteCarloStandard : public MonteCarlo
{
public:
    MonteCarloStandard(WaveFunction* wave, int number_particles, int dimension, double charge, int rank, double step_length);
    void sample(int max_variations, int number_cycles, double *cumulative_e, double *cumulative_e2, double *all_energies);
    double local_energy(double **r, double alpha, double wfold);
private:
    int number_particles;
    int dimension;
    double charge;
    int rank;
    double step_length;
};

#endif // MONTECARLOSTANDARD_H
