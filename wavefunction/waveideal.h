#ifndef WAVEFUNCTIONIDEAL_H
#define WAVEFUNCTIONIDEAL_H

#include "wavefunction.h"
#include <armadillo>
using namespace arma;
#include <iostream>

class WaveIdeal : public WaveFunction
{
public:
    WaveIdeal(Config *config);
    double evaluate(vec2 r[]);
    double laplace(vec2 r[]);
    WaveFunction* clone() {
        std::cerr << "Clone not implemented for Ideal" << std::endl;
        exit(978);
        return 0;
    }
private:
};

#endif // WAVEFUNCTIONIDEAL_H
