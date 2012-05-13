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
#include "../inih/cpp/INIReader.h"

WaveFunction::WaveFunction(Config *config) :
    config(config),
    nParticles(config->nParticles()),
    nDimensions(config->nDimensions()),
    previousEvaluation(0),
    currentEvaluation(0),
    useAnalyticalLaplace(false),
    useAnalyticalGradient(false)
{
    // allocate matrices which contain the position of the particles
    // the function matrix is defined in the progam library
    rPlus = new vec2[ nParticles];
    rMinus = new vec2[ nParticles];

    rNew = new vec2[nParticles];
    rOld = new vec2[nParticles];
}

void WaveFunction::loadConfiguration(INIReader *settings) {
    useAnalyticalLaplace = settings->GetBoolean("Wave", "useAnalyticalLaplace", false);
    useAnalyticalGradient = settings->GetBoolean("Wave", "useAnalyticalGradient", false);
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

void WaveFunction::gradientNumerical(vec2 r[], vec &rGradient)
{
    double waveEvaluation = evaluate(r);
    for (int i = 0; i < nParticles; i++) {
            rPlus[i] = rMinus[i] = r[i];
            rPlus[i] = rMinus[i] = r[i];
    }
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            rPlus[i][j] = r[i][j]+h;
            rMinus[i][j] = r[i][j]-h;
            double wfminus = evaluate(rMinus);
            double wfplus  = evaluate(rPlus);
            rGradient[i * nDimensions + j] = (wfplus - wfminus)/(2*h);
            rPlus[i][j] = r[i][j];
            rMinus[i][j] = r[i][j];
        }
    }
    rGradient = rGradient / waveEvaluation;
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

///*!
//  Tells the system about which particle was moved last
//*/
//void WaveFunction::setPreviousMovedParticle(int particleNumber)
//{
//    assert(particleNumber < nParticles);
//    previousMovedParticle = particleNumber;
//}

/*!
  This function does a ratio by evaluating the system in its current
  position and dividing by the previous evaluation. This is strongly
  optimized in the subclasses and this function is only a convenience
  to make the code run gracefully also without optimization.

  @param r An array of vectors denoting all the particle positions.
*/
double WaveFunction::ratio(vec2 &particlePosition, int particleNumber)
{
    rNew[particleNumber] = particlePosition;
    currentEvaluation = evaluate(rNew);
    double ratio = currentEvaluation / previousEvaluation;
    return ratio;
}

void WaveFunction::acceptEvaluation(int movedParticle) {
    (void)movedParticle;
    previousEvaluation = currentEvaluation;
    for(int i = 0; i < nParticles; i++) {
        rOld[i] = rNew[i];
    }
}

void WaveFunction::refuseEvaluation() {
    currentEvaluation = previousEvaluation;
    for(int i = 0; i < nParticles; i++) {
        rNew[i] = rOld[i];
    }
}

void WaveFunction::initialize(vec2 r[]) {
    previousEvaluation = evaluate(r);
    for(int i = 0; i < nParticles; i++) {
        rNew[i] = r[i];
        rOld[i] = r[i];
    }
}

WaveFunction::~WaveFunction()
{
    delete [] rPlus;
    delete [] rMinus;
}
