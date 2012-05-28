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
    void setAliveOld(bool arg) {
        m_aliveOld = arg;
    }

    bool aliveOld() {
        return m_aliveOld;
    }
    void setTimeStep(double arg) {
        timeStep = arg;
    }

    double getTimeStep() {
        return timeStep;
    }

    void progressToNextStep();

    void advance(double trialEnergy);
private:
    DiffusionWalker** otherWalkers;
    int nWalkersMax;

    double diffConstant;
    double timeStep;
    bool m_aliveNew;
    bool m_aliveOld;
};

#endif // DIFFUSIONWALKER_H
