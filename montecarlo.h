#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "wavefunction.h"

class MonteCarlo
{
public:
    MonteCarlo(WaveFunction *wave);

    void sample(int max_variations, int number_cycles, double *cumulative_e, double *cumulative_e2, double *all_energies) {};
protected:
    WaveFunction *wave;
};

#endif // MONTECARLO_H
