#ifndef WAVESIMPLE_H
#define WAVESIMPLE_H

#include "wavefunction.h"

class WaveSimple : public WaveFunction
{
public:
    WaveSimple(int m_nParticles, int m_nDimensions);
    double wave(double **r);
    double laplace(double **r);
    void setUseAnalyticalLaplace(bool val);
private:
    bool useAnalytical;
    double **rPlus;
    double **rMinus;
};

#endif // WAVESIMPLE_H
