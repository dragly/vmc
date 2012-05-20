#ifndef DIFFUSIONMONTECARLO_H
#define DIFFUSIONMONTECARLO_H

#include "montecarlo.h"
#include "../walker/diffusionwalker.h"

class DiffusionMonteCarlo : public MonteCarlo
{
public:
    DiffusionMonteCarlo(Config *config);

    void sample(int nCycles);
//    vec2 **rOld;
//    vec2 **rNew;
//    WaveFunction **waves;

    DiffusionWalker **walkers;

//    bool *aliveOld;
//    bool *aliveMid;
//    bool *aliveNew;

    int nWalkersMax;
    int nWalkersIdeal;
    int correlationStep;
    int nWalkersAlive;

    // TODO replace these names with their proper meanings
    double tau;
    double eta;
};

#endif // DIFFUSIONMONTECARLO_H
