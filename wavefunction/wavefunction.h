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
    virtual void gradient(vec2 r[], int movedParticle, vec &rGradient) {
        (void)movedParticle;
        gradientNumerical(r, rGradient);
    }
    virtual double laplace(vec2 r[], int movedParticle) {
        (void)movedParticle;
        return laplaceNumerical(r);
    }
    virtual void loadConfiguration(INIReader *settings);
    virtual double laplaceNumerical(vec2 r[]);
    virtual void gradientNumerical(vec2 r[], vec &rGradient);
    virtual void setParameters(double *parameters);
    static WaveFunction* fromName(string waveClass, Config *config);
//    virtual void setPreviousMovedParticle(int particleNumber);
    virtual double ratio(vec2 &particlePosition, int particleNumber);
    virtual void acceptEvaluation(int movedParticle);
    virtual void refuseEvaluation();
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
