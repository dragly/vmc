#include <math.h>

#include "hamiltonianstandard.h"
#include "wavefunction.h"
#include "matrix.h"
#include "utils.h"

HamiltonianStandard::HamiltonianStandard(int number_particles, int dimension, double charge) :
    number_particles(number_particles),
    dimension(dimension),
    charge(charge)
{
}

double HamiltonianStandard::energy(WaveFunction *wave, double **r, double wfold)
{
    int i, j , k;
    double e_local, wfminus, wfplus, e_kinetic, e_potential, r_12,
            r_single_particle;
    double **r_plus, **r_minus;

    // allocate matrices which contain the position of the particles
    // the function matrix is defined in the progam library
    r_plus = (double **) matrix( number_particles, dimension, sizeof(double));
    r_minus = (double **) matrix( number_particles, dimension, sizeof(double));
    for (i = 0; i < number_particles; i++) {
        for ( j=0; j < dimension; j++) {
            r_plus[i][j] = r_minus[i][j] = r[i][j];
        }
    }
    // compute the kinetic energy
    e_kinetic = 0;
    for (i = 0; i < number_particles; i++) {
        for (j = 0; j < dimension; j++) {
            r_plus[i][j] = r[i][j]+h;
            r_minus[i][j] = r[i][j]-h;
            wfminus = wave->wave(r_minus);
            wfplus  = wave->wave(r_plus);
            e_kinetic -= (wfminus+wfplus-2*wfold);
            r_plus[i][j] = r[i][j];
            r_minus[i][j] = r[i][j];
        }
    }
    // include electron mass and hbar squared and divide by wave function
    e_kinetic = 0.5*h2*e_kinetic/wfold;
    // compute the potential energy
    e_potential = 0;
    // contribution from electron-proton potential
    for (i = 0; i < number_particles; i++) {
        r_single_particle = 0;
        for (j = 0; j < dimension; j++) {
            r_single_particle += r[i][j]*r[i][j];
        }
        e_potential -= charge/sqrt(r_single_particle);
    }
    // contribution from electron-electron potential
    for (i = 0; i < number_particles-1; i++) {
        for (j = i+1; j < number_particles; j++) {
            r_12 = 0;
            for (k = 0; k < dimension; k++) {
                r_12 += (r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
            }
            e_potential += 1/sqrt(r_12);
        }
    }
    free_matrix((void **) r_plus); // free memory
    free_matrix((void **) r_minus);
    e_local = e_potential+e_kinetic;
    return e_local;
}
