#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <string>

using namespace std;

class INIReader;
class Config;

class WaveFunction
{
public:
    WaveFunction(int nParticles, int nDimensions);
    virtual double wave(double **r) = 0;
    virtual double gradient(double **r) = 0;
    virtual double laplace(double **r) = 0;
    virtual void loadConfiguration(INIReader *settings) {
        (void)settings;
    }
    double laplaceNumerical(double **r);
    void setParameters(double alpha, double beta);
    static WaveFunction* fromName(string waveClass, Config *config);
    ~WaveFunction();
protected:
    int m_nParticles;
    int m_nDimensions;
    double alpha;
    double beta;
    double **rPlus;
    double **rMinus;
};

#endif // WAVEFUNCTION_H
