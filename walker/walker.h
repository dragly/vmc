#ifndef WALKER_H
#define WALKER_H

#include <armadillo>
#include "../config.h"

using namespace arma;

class Walker
{
//    friend class Walker;
public:
    Walker(Config *config);

    virtual void advance() {}
    virtual void progressToNextStep();
    virtual void initialize(vec2 *positions);

    virtual void copyFromOther(Walker *otherWalker);

    vec2 *positionsNew() {
        return rNew;
    }

protected:
    vec2 *rNew;
    vec2 *rOld;

    vec quantumForceOld;
    vec quantumForceNew;


    double localEnergyNew;
    double localEnergyOld;

    int nParticles;
    int nDimensions;
    WaveFunction *wave;
    Hamiltonian *hamiltonian;
};

#endif // WALKER_H
