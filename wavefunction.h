#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

class WaveFunction
{
public:
    WaveFunction();
    virtual double wave(double **r) = 0;
    void setParameters(double alpha, double beta);
protected:
    double alpha;
    double beta;
};

#endif // WAVEFUNCTION_H
