#ifndef WAVESTANDARD_H
#define WAVESTANDARD_H
#include <armadillo>
using namespace arma;

#include "wavefunction.h"

class WaveStandard : public WaveFunction
{
public:
    WaveStandard(int number_particles, int dimension);
    double wave(vec2 *r);
private:
    int number_particles;
    int dimension;
};

#endif // WAVESTANDARD_H
