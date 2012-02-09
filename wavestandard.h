#ifndef WAVESTANDARD_H
#define WAVESTANDARD_H

#include "wavefunction.h"

class WaveStandard : public WaveFunction
{
public:
    WaveStandard(int number_particles, int dimension);
    double wave(double **r);
private:
    int number_particles;
    int dimension;
};

#endif // WAVESTANDARD_H
