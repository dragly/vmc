#include <math.h>

#include "wavestandard.h"

WaveStandard::WaveStandard(int number_particles, int dimension)
{
    this->number_particles = number_particles;
    this->dimension = dimension;
}

double WaveStandard::wave(double **r)
{
    int i, j;
    double wf, argument, r_single_particle;

    argument = wf = 0;
    for (i = 0; i < number_particles; i++) {
        r_single_particle = 0;
        for (j = 0; j < dimension; j++) {
            r_single_particle  += r[i][j]*r[i][j];
        }
        argument += sqrt(r_single_particle);
    }
    wf = exp(-argument*alpha) ;
    return wf;
}
