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
    virtual void gradient(const vec2 &r, int particleNumber, vec2 &rGradient) {
        gradientNumerical(r, particleNumber, rGradient);
    }
    virtual double laplace(vec2 r[]) {
        return laplaceNumerical(r);
    }
    virtual void loadConfiguration(INIReader *settings);
    virtual double laplaceNumerical(vec2 r[]);
    virtual void gradientNumerical(const vec2 &r, int particleNumber, vec2 &rGradient);
    virtual void setParameters(double *parameters);
    static WaveFunction* fromName(string waveClass, Config *config);
//    virtual void setPreviousMovedParticle(int particleNumber);
    virtual double ratio(vec2 &rParticle, int particleNumber);
    virtual void acceptEvaluation(int movedParticle);
    virtual void initialize(vec2 r[]);
    void setUseAnalyticalLaplace(bool val) {
        useAnalyticalLaplace = val;
    }

    void setUseAnalyticalGradient(bool val) {
        useAnalyticalGradient = val;
    }

    ~WaveFunction();
protected:
    Config *config;
    int previousMovedParticle;
    int nParticles;
    int nDimensions;
    double *parameters;
    vec2 *rPlus;
    vec2 *rMinus;

    vec2 *rNew;
    vec2 *rOld;
    double previousEvaluation;
    double currentEvaluation;

    bool useAnalyticalLaplace;
    bool useAnalyticalGradient;
};

#endif // WAVEFUNCTION_H
