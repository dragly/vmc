#include <string>
#include <stdio.h>

#include <armadillo>

#include "wavefunction.h"

#include "../matrix.h"
#include "../utils.h"

#include "wavesimple.h"
#include "waveideal.h"
#include "waveslater.h"
#include "../config.h"

WaveFunction::WaveFunction(Config *config) :
    m_config(config),
    m_nParticles(config->nParticles()),
    m_nDimensions(config->nDimensions())
{
    // allocate matrices which contain the position of the particles
    // the function matrix is defined in the progam library
    rPlus = new vec2[ m_nParticles];
    rMinus = new vec2[ m_nParticles];
}

double WaveFunction::laplaceNumerical(vec2 r[])
{
    double eKinetic = 0;
    double wfold = wave(r);
    for (int i = 0; i < m_nParticles; i++) {
        for (int j=0; j < m_nDimensions; j++) {
            rPlus[i][j] = rMinus[i][j] = r[i][j];
        }
    }
    for (int i = 0; i < m_nParticles; i++) {
        for (int j = 0; j < m_nDimensions; j++) {
            rPlus[i][j] = r[i][j]+h;
            rMinus[i][j] = r[i][j]-h;
            double wfminus = wave(rMinus);
            double wfplus  = wave(rPlus);
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
    for (int j = 0; j < m_nDimensions; j++) {
        for (int i = 0; i < m_nParticles; i++) {
            rPlus[i][j] = rMinus[i][j] = r[i][j];
            rPlus[i][j] = rMinus[i][j] = r[i][j];
        }
    }
    for (int j = 0; j < m_nDimensions; j++) {
        rGradient[j] = 0;
        for (int i = 0; i < m_nParticles; i++) {
            rPlus[i][j] = r[i][j]+h;
            rMinus[i][j] = r[i][j]-h;
            double wfminus = wave(rMinus);
            double wfplus  = wave(rPlus);
            rGradient[j] += (wfplus - wfminus)/(2*h);
            rPlus[i][j] = r[i][j];
            rMinus[i][j] = r[i][j];
        }
    }
}

void WaveFunction::setParameters(double* parameters)
{
    this->m_parameters = parameters;
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


WaveFunction::~WaveFunction()
{
    delete [] rPlus;
    delete [] rMinus;
}
