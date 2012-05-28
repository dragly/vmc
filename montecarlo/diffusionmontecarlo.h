#ifndef DIFFUSIONMONTECARLO_H
#define DIFFUSIONMONTECARLO_H

#include "montecarlo.h"
#include "../walker/diffusionwalker.h"

class DiffusionMonteCarlo : public MonteCarlo
{
public:
    DiffusionMonteCarlo(Config *config);

    void sample(int nSamplesLocal);
    void sample();
//    vec2 **rOld;
//    vec2 **rNew;
//    WaveFunction **waves;

    DiffusionWalker **walkers;

//    bool *aliveOld;
//    bool *aliveMid;
//    bool *aliveNew;

    int correlationStep;
    int nWalkersMax;
    int nWalkersIdeal;
    int nWalkersAlive;
    int nSamples;
    int nThermalizationCycles;
    MonteCarlo *initialMonteCarlo;
    double parameters[2];
    void loadConfiguration(INIParser *settings);
    void limitAliveWalkers();
};

#endif // DIFFUSIONMONTECARLO_H
