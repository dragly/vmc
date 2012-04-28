#ifndef WAVESIMPLE_H
#define WAVESIMPLE_H

#include <armadillo>

#include "wavefunction.h"

using namespace std;
using namespace arma;

class WaveSimple : public WaveFunction
{
public:
    WaveSimple(Config *config);
    double wave(vec2 r[]);
    double laplace(vec2 r[]);
    void setUseAnalyticalLaplace(bool val);
private:
    bool useAnalytical;
    vec2 rPlus[];
    vec2 rMinus[];
};

#endif // WAVESIMPLE_H
