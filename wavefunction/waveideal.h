#ifndef WAVEFUNCTIONIDEAL_H
#define WAVEFUNCTIONIDEAL_H
#include <armadillo>
using namespace arma;

#include "wavefunction.h"

class WaveIdeal : public WaveFunction
{
public:
    WaveIdeal(Config *config);
    double wave(vec2 *r);
    double laplace(vec2 *r);
    void setUseAnalyticalLaplace(bool val){
        useAnalytical = val;
    }
private:
    bool useAnalytical;
};

#endif // WAVEFUNCTIONIDEAL_H
