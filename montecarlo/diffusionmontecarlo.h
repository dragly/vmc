#ifndef DIFFUSIONMONTECARLO_H
#define DIFFUSIONMONTECARLO_H

#include "montecarlo.h"
#include "../walker/diffusionwalker.h"

/*!
  * \brief Implements the diffusion Monte Carlo algorithm
  */
class DiffusionMonteCarlo : public MonteCarlo
{
public:
    DiffusionMonteCarlo(Config *config);

    void sample(int nSamplesLocal);
    void sample();
    void setTimeStep(double arg) {
        timeStep = arg;
        updateWalkerParameters();
    }
    void loadConfiguration(INIParser *settings);

//    vec2 **rOld;
//    vec2 **rNew;
//    WaveFunction **waves;
private:
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
    double timeStep;
    double parameters[2];
    void updateWalkerParameters();
};

#endif // DIFFUSIONMONTECARLO_H
