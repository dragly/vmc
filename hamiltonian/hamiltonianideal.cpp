#include "hamiltonian/hamiltonianideal.h"
#include "matrix.h"
#include "utils.h"
#include "wavefunction.h"
#include <math.h>

HamiltonianIdeal::HamiltonianIdeal(int number_particles, int dimension, double charge) :
    number_particles(number_particles),
    dimension(dimension),
    charge(charge)
{
}

double HamiltonianIdeal::energy(WaveFunction *wave, double **r)
{
    double eLocal, eKinetic, ePotential,
            rSingleParticle;
    // compute the kinetic energy
    // TODO: Use the exact derivative
    // TODO: Create a derivative-finder function that uses interpolation to approximate the derivative
    eKinetic = -0.5*wave->laplace(r);
    //    e_kinetic = 0.5*e_kinetic/h2;
    // compute the potential energy
    ePotential = 0;
    // contribution from harmonic oscillator potential
    for (int i = 0; i < number_particles; i++) {
        rSingleParticle = 0;
        for (int j = 0; j < dimension; j++) {
            rSingleParticle += r[i][j]*r[i][j];
        }
        ePotential += 0.5 * rSingleParticle;
    }
    //     contribution from electron-electron potential
    for (int i = 0; i < number_particles-1; i++) {
        for (int j = i+1; j < number_particles; j++) {
            double r_12 = 0;
            for (int k = 0; k < dimension; k++) {
                r_12 += (r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
            }
            ePotential += 1/sqrt(r_12);
        }
    }
    eLocal = ePotential+eKinetic;
    return eLocal;
}

double HamiltonianIdeal::kineticEnergy(WaveFunction *wave, double **r, double wfold)
{
    double **r_plus, **r_minus;

    // allocate matrices which contain the position of the particles
    // the function matrix is defined in the progam library
    r_plus = (double **) matrix( number_particles, dimension, sizeof(double));
    r_minus = (double **) matrix( number_particles, dimension, sizeof(double));
    for (int i = 0; i < number_particles; i++) {
        for (int j=0; j < dimension; j++) {
            r_plus[i][j] = r_minus[i][j] = r[i][j];
        }
    }
    double e_kinetic = 0;
    for (int i = 0; i < number_particles; i++) {
        for (int j = 0; j < dimension; j++) {
            r_plus[i][j] = r[i][j]+h;
            r_minus[i][j] = r[i][j]-h;
            double wfminus = wave->wave(r_minus);
            double wfplus  = wave->wave(r_plus);
            e_kinetic -= (wfminus+wfplus-2*wfold);
            r_plus[i][j] = r[i][j];
            r_minus[i][j] = r[i][j];
        }
    }
    // include electron mass and hbar squared and divide by wave function
    e_kinetic = 0.5*h2*e_kinetic/wfold;

    free_matrix((void **) r_plus); // free memory
    free_matrix((void **) r_minus);
    return e_kinetic;
}
