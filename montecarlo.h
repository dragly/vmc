#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "wavefunction.h"
#include "hamiltonian/hamiltonian.h"

class MonteCarlo
{
public:
    MonteCarlo(WaveFunction *wave, Hamiltonian *hamiltonian);

    virtual void sample(int number_cycles, double *energies) = 0;
protected:
    WaveFunction *wave;
    Hamiltonian *hamiltonian;
};

#endif // MONTECARLO_H
