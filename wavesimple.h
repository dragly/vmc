#ifndef WAVESIMPLE_H
#define WAVESIMPLE_H

#include "wavefunction.h"

class WaveSimple : public WaveFunction
{
public:
    WaveSimple(int number_particles, int dimension);
    double wave(double **r);
private:
    int number_particles;
    int dimension;
};

#endif // WAVESIMPLE_H
