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
    charge(charge)
{
}

double HamiltonianSimple::energy(WaveFunction *wave, double **r, double alpha, double wfold)
{
    int i, j;
    double e_local, e_kinetic, e_potential,
            r_single_particle;
    // compute the kinetic energy
    // TODO: Use the exact derivative
    // TODO: Create a derivative-finder function that uses interpolation to approximate the derivative
    e_kinetic = kineticEnergy(wave, r, alpha, wfold);
    //    e_kinetic = 0.5*e_kinetic/h2;
    // compute the potential energy
    e_potential = 0;
    // contribution from electron-proton potential
    for (i = 0; i < number_particles; i++) {
        r_single_particle = 0;
        for (j = 0; j < dimension; j++) {
            r_single_particle += r[i][j]*r[i][j];
        }
        e_potential += 0.5 * r_single_particle;
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
    e_local = e_potential+e_kinetic;
    return e_local;
}

double HamiltonianSimple::kineticEnergy(WaveFunction *wave, double **r, double alpha, double wfold)
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
            double wfminus = wave->wave(r_minus, alpha);
            double wfplus  = wave->wave(r_plus, alpha);
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

/*
** The function
**           splint()
** takes xa[0,..,n - 1] and y[0,..,n - 1] which tabulates a function
** (with the xa[i]'s in order) and given ya[0,..,n - 1], which is the
** output from function spline() and with given value of x returns a
** cubic--spline interpolation value y.
**
** The function is borrowed from the lib.cpp library from the FYS3150 course
** by Morten H. Jensen
*/

void HamiltonianSimple::splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
    int       klo,khi,k;
    double    ha,b,a;

    klo = 0;
    khi = n - 1;
    while((khi - klo) > 1) {   // binary search
        k = (khi + klo) >> 1;
        if(xa[k] > x)   khi = k;
        else            klo = k;
    }
    ha = xa[khi] - xa[klo];
    if(fabs(ha) < ZERO) {
        printf("\n\n Error in function splint(): ");
        printf("\n The difference h = %4.1E -- too small\n",h);
        exit(1);
    }
    a  = (xa[khi] - x)/ha;
    b  = (x - xa[klo])/ha;
    *y =   a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo]
                                        + (b * b * b - b) * y2a[khi]) * (ha * ha)/6.0;

} // End: function splint()
