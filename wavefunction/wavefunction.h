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
    virtual double evaluate(vec2 r[]) = 0;
    virtual void gradient(vec2 r[], vec2 &rGradient) {
        gradientNumerical(r, rGradient);
    }
    virtual double laplace(vec2 r[]) {
        return laplaceNumerical(r);
    }
    virtual void loadConfiguration(INIReader *settings) {
        (void)settings;
    }
    virtual double laplaceNumerical(vec2 r[]);
    virtual void gradientNumerical(vec2 r[], vec2 &rGradient);
    virtual void setParameters(double *parameters);
    static WaveFunction* fromName(string waveClass, Config *config);
    virtual void setPreviousMovedParticle(int particleNumber);
    virtual double ratio(vec2 rNew[]);
    virtual void acceptEvaluation();
    virtual void init(vec2 r[]);

    ~WaveFunction();
protected:
    Config *config;
    int previousMovedParticle;
    int nParticles;
    int nDimensions;
    double *parameters;
    vec2 *rPlus;
    vec2 *rMinus;
    double previousEvaluation;
    double currentEvaluation;
};

#endif // WAVEFUNCTION_H
