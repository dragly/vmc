#ifndef DIFFUSIONMONTECARLO_H
#define DIFFUSIONMONTECARLO_H

#include "montecarlo.h"

class DiffusionMonteCarlo : public MonteCarlo
{
public:
    DiffusionMonteCarlo(Config *config);

    void sample(int numberCycles);
};

#endif // DIFFUSIONMONTECARLO_H
