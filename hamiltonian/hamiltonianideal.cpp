#include "hamiltonianideal.h"
#include "../matrix.h"
#include "../utils.h"
#include "../wavefunction/wavefunction.h"
#include "../config.h"
#include <math.h>

HamiltonianIdeal::HamiltonianIdeal(Config *config) :
    Hamiltonian(config),
    omega(config->omega())
{
}
int ntries = 0;

double HamiltonianIdeal::kineticEnergy(WaveFunction *wave, vec2 r[])
{
    return -0.5*wave->laplace(r);
}

double HamiltonianIdeal::potentialEnergy(WaveFunction *wave, vec2 r[])
{
    double ePotential = 0;
    ePotential += externalPotentialEnergy(wave, r);
    if(m_interactionEnabled) {
        ePotential += interactionPotentialEnergy(wave, r);
    }
    return ePotential;
}

double HamiltonianIdeal::externalPotentialEnergy(WaveFunction *wave, vec2 r[])
{
    (void)wave;
    double externalEnergy = 0;
    double rSingleParticle = 0;
    // contribution from harmonic oscillator potential
    for (int i = 0; i < m_nParticles; i++) {
        rSingleParticle = 0;
        for (int j = 0; j < m_nDimensions; j++) {
            rSingleParticle += r[i][j]*r[i][j];
        }
        externalEnergy += 0.5 * omega * omega * rSingleParticle;
    }
    return externalEnergy;
}

double HamiltonianIdeal::interactionPotentialEnergy(WaveFunction *wave, vec2 r[])
{
    (void)wave;
    double interactionEnergy = 0;
    //     contribution from electron-electron potential
    // TODO Optimization - store an array of potential energies between particles and only update necessary parts (reduces use of sqrt)
    for (int i = 0; i < m_nParticles-1; i++) {
        for (int j = i+1; j < m_nParticles; j++) {
            double r_12 = 0;
            for (int k = 0; k < m_nDimensions; k++) {
                r_12 += (r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
            }
            interactionEnergy += 1/sqrt(r_12);
        }
    }
    return interactionEnergy;
}
