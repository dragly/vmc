#include <string>
#include <stdio.h>

#include <armadillo>

#include <assert.h>

#include "wavefunction.h"

#include "../matrix.h"
#include "../utils.h"

#include "wavesimple.h"
#include "waveideal.h"
#include "waveslater.h"
#include "../config.h"

WaveFunction::WaveFunction(Config *config) :
    config(config),
    nParticles(config->nParticles()),
    nDimensions(config->nDimensions()),
    previousEvaluation(0),
    currentEvaluation(0)
{
    // allocate matrices which contain the position of the particles
    // the function matrix is defined in the progam library
    rPlus = new vec2[ nParticles];
    rMinus = new vec2[ nParticles];
}

double WaveFunction::laplaceNumerical(vec2 r[])
{
    double eKinetic = 0;
    double wfold = evaluate(r);
    for (int i = 0; i < nParticles; i++) {
        for (int j=0; j < nDimensions; j++) {
            rPlus[i][j] = rMinus[i][j] = r[i][j];
        }
    }
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            rPlus[i][j] = r[i][j]+h;
            rMinus[i][j] = r[i][j]-h;
            double wfminus = evaluate(rMinus);
            double wfplus  = evaluate(rPlus);
            eKinetic += h2*(wfminus+wfplus-2*wfold);
            rPlus[i][j] = r[i][j];
            rMinus[i][j] = r[i][j];
        }
    }
    eKinetic /= wfold;

    return eKinetic;
}

void WaveFunction::gradientNumerical(vec2 r[], vec2 &rGradient)
{
    for (int j = 0; j < nDimensions; j++) {
        for (int i = 0; i < nParticles; i++) {
            rPlus[i][j] = rMinus[i][j] = r[i][j];
            rPlus[i][j] = rMinus[i][j] = r[i][j];
        }
    }
    for (int j = 0; j < nDimensions; j++) {
        rGradient[j] = 0;
        for (int i = 0; i < nParticles; i++) {
            rPlus[i][j] = r[i][j]+h;
            rMinus[i][j] = r[i][j]-h;
            double wfminus = evaluate(rMinus);
            double wfplus  = evaluate(rPlus);
            rGradient[j] += (wfplus - wfminus)/(2*h);
            rPlus[i][j] = r[i][j];
            rMinus[i][j] = r[i][j];
        }
    }
}

void WaveFunction::setParameters(double* parameters)
{
    this->parameters = parameters;
}

/*!
  Creates and returns a new wave function based on the class specified by the parameter string.
  */
WaveFunction* WaveFunction::fromName(std::string waveClass, Config* config) {
    if(waveClass == "WaveSimple") {
        WaveSimple *waveSimple = new WaveSimple(config);
        return waveSimple;
    } else if(waveClass == "WaveIdeal") {
        WaveIdeal *waveIdeal = new WaveIdeal(config);
        return waveIdeal;
    } else if(waveClass == "WaveSlater") {
        WaveSlater *waveIdeal = new WaveSlater(config);
        return waveIdeal;
    } else {
        return 0;
    }
}

/*!
  Tells the system about which particle was moved last
*/
void WaveFunction::setPreviousMovedParticle(int particleNumber)
{
    assert(particleNumber < nParticles);
    previousMovedParticle = particleNumber;
}

/*!
  This function does a ratio by evaluating the system in its current
  position and dividing by the previous evaluation. This is strongly
  optimized in the subclasses and this function is only a convenience
  to make the code run gracefully also without optimization.

  @param r An array of vectors denoting all the particle positions.
*/
double WaveFunction::ratio(vec2 rNew[])
{
    currentEvaluation = evaluate(rNew);
    double ratio = currentEvaluation / previousEvaluation;
    return ratio;
}

void WaveFunction::acceptEvaluation() {
    previousEvaluation = currentEvaluation;
}

void WaveFunction::init(vec2 r[]) {
    previousEvaluation = evaluate(r);
}

WaveFunction::~WaveFunction()
{
    delete [] rPlus;
    delete [] rMinus;
}
