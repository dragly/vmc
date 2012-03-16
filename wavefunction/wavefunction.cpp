#include <string>
#include <stdio.h>

#include <armadillo>

#include "wavefunction.h"

#include "../matrix.h"
#include "../utils.h"

#include "wavesimple.h"
#include "waveideal.h"
#include "../config.h"

WaveFunction::WaveFunction(int nParticles, int nDimensions) :
    m_nParticles(nParticles),
    m_nDimensions(nDimensions)
{
    // allocate matrices which contain the position of the particles
    // the function matrix is defined in the progam library
    rPlus = new vec2[ nParticles];
    rMinus = new vec2[ nParticles];
}

double WaveFunction::laplaceNumerical(vec2 *r)
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

void WaveFunction::gradientNumerical(vec2 *r, double *rGradient)
{
    for (int i = 0; i < m_nParticles; i++) {
        for (int j=0; j < m_nDimensions; j++) {
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

void WaveFunction::setParameters(double alpha, double beta)
{
    this->alpha = alpha;
    this->beta = beta;
}

/*!
  Creates and returns a new wave function based on the class specified by the parameter string.
  */
WaveFunction* WaveFunction::fromName(std::string waveClass, Config* config) {
    if(waveClass == "WaveSimple") {
        WaveSimple *waveSimple = new WaveSimple(config->nParticles(), config->nDimensions());
        return waveSimple;
    } else if(waveClass == "WaveIdeal") {
        WaveIdeal *waveIdeal = new WaveIdeal(config->nParticles(), config->nDimensions());
        return waveIdeal;
    } else {
        return 0;
    }
}


WaveFunction::~WaveFunction()
{
    free_matrix((void **) rPlus); // free memory
    free_matrix((void **) rMinus);
}
