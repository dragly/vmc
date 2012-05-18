#ifndef DIFFUSIONMONTECARLO_H
#define DIFFUSIONMONTECARLO_H

#include "montecarlo.h"

class DiffusionMonteCarlo : public MonteCarlo
{
public:
    DiffusionMonteCarlo(Config *config);

    void sample(int nCycles);
    vec2 **rOld;
    vec2 **rNew;
    WaveFunction **waves;

    int nTotalWalkers;
    int correlationStep;

    // TODO replace these names with their proper meanings
    double tau;
    double eta;
};

#endif // DIFFUSIONMONTECARLO_H
