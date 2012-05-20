#include "walker.h"

Walker::Walker(Config *config) :
    nParticles(config->nParticles()),
    nDimensions(config->nDimensions()),
    wave(config->wave()),
    hamiltonian(config->hamiltonian())
{
    rNew = new vec2[nParticles];
    rOld = new vec2[nParticles];

    quantumForceNew = zeros<vec>(nDimensions * nParticles);
    quantumForceOld = zeros<vec>(nDimensions * nParticles);

}

void Walker::progressToNextStep()
{
    for(int i = 0; i < nParticles; i++) {
        rOld[i] = rNew[i];
    }
    localEnergyOld = localEnergyNew;
    quantumForceOld = quantumForceNew;
}

void Walker::initialize(vec2 *positions)
{
    for(int i = 0; i < nParticles; i++) {
        rNew[i] = positions[i];
        rOld[i] = positions[i];
    }
    wave->initialize(positions);
    wave->gradient(positions, 0, quantumForceNew);
    quantumForceOld = quantumForceNew;
    localEnergyOld = hamiltonian->energy(wave, rNew);
    localEnergyNew = localEnergyOld;
}

void Walker::copyFromOther(Walker *otherWalker) {
    for(int i = 0; i < nParticles; i++) {
        rNew[i] = otherWalker->rNew[i];
        rOld[i] = otherWalker->rOld[i];
    }
    quantumForceNew = otherWalker->quantumForceNew;
    quantumForceOld = otherWalker->quantumForceOld;
    wave->initialize(otherWalker->rNew);
    localEnergyNew = otherWalker->localEnergyNew;
    localEnergyOld = otherWalker->localEnergyOld;
}
