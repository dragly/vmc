#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "hamiltoniansimple.h"
#include "../wavefunction/wavefunction.h"
#include "../matrix.h"
#include "../utils.h"

HamiltonianSimple::HamiltonianSimple(Config *config) :
    Hamiltonian(config)
{
}

double HamiltonianSimple::potentialEnergy(WaveFunction *wave, vec2 r[])
{
    (void)wave;
    double rSingleParticle;
    double ePotential = 0;
    // contribution from electron-proton potential
    for (int i = 0; i < m_nParticles; i++) {
        rSingleParticle = 0;
        for (int j = 0; j < m_nDimensions; j++) {
            rSingleParticle += r[i][j]*r[i][j];
        }
        ePotential += 0.5 * rSingleParticle;
    }
    return ePotential;
}

double HamiltonianSimple::kineticEnergy(WaveFunction *wave, vec2 r[])
{
    return -0.5*wave->laplace(r);
}
