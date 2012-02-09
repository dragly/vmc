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

double HamiltonianIdeal::energy(WaveFunction *wave, double **r, double wfold)
{
    double e_local, e_kinetic, e_potential,
            r_single_particle;
    // compute the kinetic energy
    // TODO: Use the exact derivative
    // TODO: Create a derivative-finder function that uses interpolation to approximate the derivative
    e_kinetic = kineticEnergy(wave, r, wfold);
    //    e_kinetic = 0.5*e_kinetic/h2;
    // compute the potential energy
    e_potential = 0;
    // contribution from harmonic oscillator potential
    for (int i = 0; i < number_particles; i++) {
        r_single_particle = 0;
        for (int j = 0; j < dimension; j++) {
            r_single_particle += r[i][j]*r[i][j];
        }
        e_potential += 0.5 * r_single_particle;
    }
    //     contribution from electron-electron potential
    for (int i = 0; i < number_particles-1; i++) {
        for (int j = i+1; j < number_particles; j++) {
            double r_12 = 0;
            for (int k = 0; k < dimension; k++) {
                r_12 += (r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
            }
            e_potential += 1/sqrt(r_12);
        }
    }
    e_local = e_potential+e_kinetic;
    return e_local;
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
