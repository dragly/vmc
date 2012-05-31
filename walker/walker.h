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
    virtual void initialize(vec2 *positions);

    virtual void copyOtherWalker(Walker *otherWalker);

    vec2 *positionsNew() {
        return rNew;
    }
    double energy() {
        return m_energy;
    }

    int changeInEnergySamples() {
        return m_changeInEnergySamples;
    }

    void setParameters(double *parameters) {
        wave->setParameters(parameters);
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
    double m_energy;
    long *idum;
    int m_changeInEnergySamples;
};

#endif // WALKER_H
