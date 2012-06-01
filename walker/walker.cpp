#include "walker.h"

Walker::Walker(Config *config) :
    nParticles(config->nParticles()),
    nDimensions(config->nDimensions()),
    hamiltonian(config->hamiltonian()),
    m_energy(0),
    idum(config->idum()),
    isCopyFromOtherWalker(false)
{
    wave = config->wave()->clone();
//    wave = config->wave();

    rNew = new vec2[nParticles];
    rOld = new vec2[nParticles];

    quantumForceNew = zeros<vec>(nDimensions * nParticles);
    quantumForceOld = zeros<vec>(nDimensions * nParticles);

}

void Walker::initialize(vec2 *positions)
{
    for(int i = 0; i < nParticles; i++) {
        rNew[i] = positions[i];
        rOld[i] = positions[i];
    }
    wave->initialize(positions);
    wave->gradient(positions, quantumForceNew);
    quantumForceOld = quantumForceNew;
    localEnergyOld = hamiltonian->energy(wave, rNew);
    if(isnan(localEnergyOld)) {
        std::cout << rNew[0] << rNew[1] << std::endl;
    }
    localEnergyNew = localEnergyOld;
}

void Walker::copyOtherWalker(Walker *otherWalker) {
    isCopyFromOtherWalker = true;
    for(int i = 0; i < nParticles; i++) {
        rNew[i] = otherWalker->rNew[i];
        rOld[i] = otherWalker->rOld[i];
    }
    quantumForceNew = otherWalker->quantumForceNew;
    quantumForceOld = otherWalker->quantumForceOld;
    wave->initialize(rOld);
    localEnergyNew = otherWalker->localEnergyNew;
    localEnergyOld = otherWalker->localEnergyOld;
}
