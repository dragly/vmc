#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>
using namespace arma;

#include <string>
using namespace std;

class INIParser;
class Config;

class WaveFunction
{
public:
    WaveFunction(Config *config);
    virtual double evaluate(vec2 r[]) = 0;
    virtual void gradient(vec2 r[], vec &rGradient) {
        gradientNumerical(r, rGradient);
    }
    virtual double laplace(vec2 r[]) {
        return laplaceNumerical(r);
    }
    virtual WaveFunction* clone() = 0;
    virtual void loadConfiguration(INIParser *settings);
    virtual double laplaceNumerical(vec2 r[]);
    virtual void gradientNumerical(vec2 r[], vec &rGradient);
    virtual void setParameters(double *parameters);
    static WaveFunction* fromName(string waveClass, Config *config);
//    virtual void setPreviousMovedParticle(int particleNumber);
    virtual double ratio(vec2 &particlePosition, int movedParticle);
    virtual void acceptMove(int movedParticle);
    virtual void rejectMove();
    virtual void initialize(vec2 r[]);
    virtual void prepareGradient(vec2 &particlePosition, int movedParticle) {
        (void) movedParticle;
        (void) particlePosition;
    }

    void setUseAnalyticalLaplace(bool val) {
        useAnalyticalLaplace = val;
    }

    void setUseAnalyticalGradient(bool val) {
        useAnalyticalGradient = val;
    }

    ~WaveFunction();
    virtual void outputProperties();
protected:
    Config *config;
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
