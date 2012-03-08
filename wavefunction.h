#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

class INIReader;

class WaveFunction
{
public:
    WaveFunction(int number_particles, int dimension);
    virtual double wave(double **r) = 0;
    virtual double gradient(double **r) = 0;
    virtual double laplace(double **r) = 0;
    virtual void loadConfiguration(INIReader *settings) {}
    double laplaceNumerical(double **r);
    void setParameters(double alpha, double beta);
    ~WaveFunction();
protected:
    int nParticles;
    int dimensions;
    double alpha;
    double beta;
    double **rPlus;
    double **rMinus;
};

#endif // WAVEFUNCTION_H
