#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "hamiltoniansimple.h"
#include "wavefunction.h"
#include "matrix.h"
#include "utils.h"

HamiltonianSimple::HamiltonianSimple(int number_particles, int dimension, double charge) :
    number_particles(number_particles),
    dimension(dimension),
    charge(charge),
    useAnalyticalKineticEnergy(false)
{
}

double HamiltonianSimple::energy(WaveFunction *wave, double **r, double wfold)
{
    int i, j;
    double eLocal;
    double eKinetic;
    double ePotential;
    double rSingleParticle;
    // compute the kinetic energy
    // TODO: Use the exact derivative
    // TODO: Create a derivative-finder function that uses interpolation to approximate the derivative
    eKinetic = -0.5*wave->laplace(r);
    //    e_kinetic = 0.5*e_kinetic/h2;
    // compute the potential energy
    ePotential = 0;
    // contribution from electron-proton potential
    for (i = 0; i < number_particles; i++) {
        rSingleParticle = 0;
        for (j = 0; j < dimension; j++) {
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
