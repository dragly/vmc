#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "wavefunction.h"
#include "hamiltonian.h"

class MonteCarlo
{
public:
    MonteCarlo(WaveFunction *wave, Hamiltonian *hamiltonian);

    virtual void sample(int max_variations, int number_cycles, double *cumulative_e, double *cumulative_e2, double *all_energies) = 0;
protected:
    WaveFunction *wave;
    Hamiltonian *hamiltonian;
};

#endif // MONTECARLO_H
