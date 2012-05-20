#ifndef DIFFUSIONWALKER_H
#define DIFFUSIONWALKER_H

#include "walker.h"

class DiffusionWalker : public Walker
{
public:
    DiffusionWalker(Config *config, DiffusionWalker** otherWalkers, int nWalkersMax);
    bool aliveNew() {
        return m_aliveNew;
    }
    void setAliveNew(bool arg) {
        m_aliveNew = arg;
    }

    bool aliveOld() {
        return m_aliveOld;
    }
    void progressToNextStep();

    double energy() {
        return m_energy;
    }

    int changeInWalkersAlive() {
        return m_changeInWalkersAlive;
    }

    int changeInEnergySamples() {
        return m_changeInEnergySamples;
    }

    void advance(double trialEnergy);
private:
    bool m_aliveNew;
    bool m_aliveOld;
    DiffusionWalker** otherWalkers;
    int nWalkersMax;

    double diffConstant;
    double tau;
    long *idum;
    int m_changeInWalkersAlive;
    int m_changeInEnergySamples;
    double m_energy;
};

#endif // DIFFUSIONWALKER_H
