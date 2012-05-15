#ifndef WAVEFUNCTIONIDEAL_H
#define WAVEFUNCTIONIDEAL_H
#include <armadillo>
using namespace arma;

#include "wavefunction.h"

class WaveIdeal : public WaveFunction
{
public:
    WaveIdeal(Config *config);
    double evaluate(vec2 r[]);
    double laplace(vec2 r[], int movedParticle);
private:
};

#endif // WAVEFUNCTIONIDEAL_H
