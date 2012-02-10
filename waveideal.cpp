#include "waveideal.h"
#include <math.h>

WaveIdeal::WaveIdeal(int number_particles, int dimension) :
    WaveFunction(),
    number_particles(number_particles),
    dimension(dimension)
{
}

double WaveIdeal::wave(double **r)
{
    int i, j;
    double wf, argument, r_single_particle;

    argument = wf = 0;
    for (i = 0; i < number_particles; i++) {
        r_single_particle = 0;
        for (j = 0; j < dimension; j++) {
            r_single_particle  += r[i][j]*r[i][j];
        }
        argument += r_single_particle;
    }
    double vec[2];
    vec[0] = r[0][0]-r[1][0];
    vec[1] = r[0][1]-r[1][1];
    double r12 = sqrt(vec[0]*vec[0]+vec[1]*vec[1]);
    double a = 1;
    wf = exp(-(argument*alpha) / 2);
    double jastrowArgument = (a * r12) / (1 + beta * r12);
    wf *= exp(jastrowArgument);
    return wf;
}

double WaveIdeal::laplace(double **r)
{
    return 0;
}
