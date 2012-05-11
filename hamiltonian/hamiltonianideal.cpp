#include "hamiltonianideal.h"
#include "../matrix.h"
#include "../utils.h"
#include "../wavefunction/wavefunction.h"
#include <math.h>

HamiltonianIdeal::HamiltonianIdeal(Config *config) :
    Hamiltonian(config)
{
}

double HamiltonianIdeal::energy(WaveFunction *wave, vec2 r[])
{
    double eLocal, eKinetic, ePotential,
            rSingleParticle;
    double omega = 1;
    // compute the kinetic energy
    // TODO: Create a derivative-finder function that uses interpolation to approximate the derivative
    // TODO: Add the number of the particle that has been moved
    eKinetic = -0.5*wave->laplace(r, 0);
    //    e_kinetic = 0.5*e_kinetic/h2;
    // compute the potential energy
    ePotential = 0;
    // contribution from harmonic oscillator potential
    for (int i = 0; i < m_nParticles; i++) {
        rSingleParticle = 0;
        for (int j = 0; j < m_nDimensions; j++) {
            rSingleParticle += r[i][j]*r[i][j];
        }
        ePotential += 0.5 * omega * omega * rSingleParticle;
    }
    if(m_interactionEnabled) {
        //     contribution from electron-electron potential
        for (int i = 0; i < m_nParticles-1; i++) {
            for (int j = i+1; j < m_nParticles; j++) {
                double r_12 = 0;
                for (int k = 0; k < m_nDimensions; k++) {
                    r_12 += (r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
                }
                ePotential += 1/sqrt(r_12);
            }
        }
    }
    eLocal = ePotential+eKinetic;
    return eLocal;
}
