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

double HamiltonianSimple::energy(WaveFunction *wave, vec2 r[])
{
    int i, j;
    double eLocal;
    double eKinetic;
    double ePotential;
    double rSingleParticle;
    // compute the kinetic energy
    // TODO: Create a derivative-finder function that uses interpolation to approximate the derivative
    // TODO: Add the number of the particle that has been moved
    eKinetic = -0.5*wave->laplace(r, 0);
    //    e_kinetic = 0.5*e_kinetic/h2;
    // compute the potential energy
    ePotential = 0;
    // contribution from electron-proton potential
    for (i = 0; i < m_nParticles; i++) {
        rSingleParticle = 0;
        for (j = 0; j < m_nDimensions; j++) {
            rSingleParticle += r[i][j]*r[i][j];
        }
        ePotential += 0.5 * rSingleParticle;
    }
    // contribution from electron-electron potential
    //    for (i = 0; i < number_particles-1; i++) {
    //        for (j = i+1; j < number_particles; j++) {
    //            r_12 = 0;
    //            for (k = 0; k < dimension; k++) {
    //                r_12 += (r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
    //            }
    //            e_potential += 1/sqrt(r_12);
    //        }
    //    }
    eLocal = ePotential+eKinetic;
    return eLocal;
}
