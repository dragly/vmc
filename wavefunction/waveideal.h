#ifndef WAVEFUNCTIONIDEAL_H
#define WAVEFUNCTIONIDEAL_H

#include "wavefunction.h"

class WaveIdeal : public WaveFunction
{
public:
    WaveIdeal(int number_particles, int dimension);
    double wave(double **r);
    double laplace(double **r);
    void setUseAnalyticalLaplace(bool val){
        useAnalytical = val;
    }
private:
    bool useAnalytical;
};

#endif // WAVEFUNCTIONIDEAL_H