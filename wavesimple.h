#ifndef WAVESIMPLE_H
#define WAVESIMPLE_H

#include "wavefunction.h"

class WaveSimple : public WaveFunction
{
public:
    WaveSimple(int number_particles, int dimension);
    double wave(double **r);
    double gradient(double **r) {return 0; }
    double laplace(double **r);
    double **r_plus, **r_minus;
    void setUseAnalytical(bool val);

    ~WaveSimple();
private:
    int number_particles;
    int dimension;
    bool useAnalytical;

    double *preExp;

    int nExp;
    double aExp;
    double bExp;
    double ratioExp;
    double ratioExpInverse;
    double nExpInverse;

    double myExp(double x);
};

#endif // WAVESIMPLE_H
