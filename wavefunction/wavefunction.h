#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include <armadillo>
using namespace arma;

#include <string>

using namespace std;

class INIReader;
class Config;

class WaveFunction
{
public:
    WaveFunction(int nParticles, int nDimensions);
    virtual double wave(vec2 *r) = 0;
    virtual void gradient(vec2 *r, vec2 &rGradient) {
        gradientNumerical(r, rGradient);
    }
    virtual double laplace(vec2 *r) = 0;
    virtual void loadConfiguration(INIReader *settings) {
        (void)settings;
    }
    double laplaceNumerical(vec2 *r);
    void gradientNumerical(vec2 *r, vec2 &rGradient);
    void setParameters(double alpha, double beta);
    static WaveFunction* fromName(string waveClass, Config *config);
    ~WaveFunction();
protected:
    int m_nParticles;
    int m_nDimensions;
    double alpha;
    double beta;
    vec2 *rPlus;
    vec2 *rMinus;
};

#endif // WAVEFUNCTION_H
