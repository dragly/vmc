#ifndef WAVEFUNCTIONIDEAL_H
#define WAVEFUNCTIONIDEAL_H

#include "wavefunction.h"

class WaveIdeal : public WaveFunction
{
public:
    WaveIdeal(int number_particles, int dimension);
    double wave(double **r);
    double gradient(double **r) { return 0; }
    double laplace(double **r);
private:
    int number_particles;
    int dimension;
};

#endif // WAVEFUNCTIONIDEAL_H
