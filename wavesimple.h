#ifndef WAVESIMPLE_H
#define WAVESIMPLE_H

#include "wavefunction.h"

class WaveSimple : public WaveFunction
{
public:
    WaveSimple(int nParticles, int dimensions);
    double wave(double **r);
    double gradient(double **r) {return 0; }
    double laplace(double **r);
    void setUseAnalyticalLaplace(bool val);
private:
    bool useAnalytical;
    double **rPlus;
    double **rMinus;
};

#endif // WAVESIMPLE_H
