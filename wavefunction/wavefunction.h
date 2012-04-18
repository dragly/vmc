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
    WaveFunction(Config *config);
    virtual double wave(const vec2 r[]) = 0;
    virtual void gradient(const vec2 r[], vec2 &rGradient) {
        gradientNumerical(r, rGradient);
    }
    virtual double laplace(const vec2 r[]) {
        return laplaceNumerical(r);
    }
    virtual void loadConfiguration(INIReader *settings) {
        (void)settings;
    }
    double laplaceNumerical(const vec2 r[]);
    void gradientNumerical(const vec2 r[], vec2 &rGradient);
    void setParameters(double *m_parameters);
    static WaveFunction* fromName(string waveClass, Config *config);
    ~WaveFunction();
protected:
    Config *m_config;
    int m_nParticles;
    int m_nDimensions;
    double *m_parameters;
    vec2 *rPlus;
    vec2 *rMinus;
};

#endif // WAVEFUNCTION_H
