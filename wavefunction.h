#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

class WaveFunction
{
public:
    WaveFunction();
    virtual double wave(double **r, double alpha) = 0;
};

#endif // WAVEFUNCTION_H
