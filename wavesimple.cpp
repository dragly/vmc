#include "wavesimple.h"
#include "matrix.h"
#include "utils.h"
#include <math.h>
#include <iostream>
using namespace std;

WaveSimple::WaveSimple(int nParticles, int dimensions)  :
    WaveFunction(nParticles, dimensions),
    useAnalytical(false)
{
}

double WaveSimple::wave(double **r)
{
    int i, j;
    double wf, argument, rSingleParticle;

    argument = wf = 0;
    for (i = 0; i < m_nParticles; i++) {
        rSingleParticle = 0;
        for (j = 0; j < m_nDimensions; j++) {
            rSingleParticle  += r[i][j]*r[i][j];
        }
        argument += rSingleParticle;
    }
    wf = exp(-(argument*alpha) / 2) ;
    return wf;
}

double WaveSimple::laplace(double **r)
{
    if(useAnalytical) {
        double eKinetic = 0;
        double rSingleParticle;
        for (int i = 0; i < m_nParticles; i++) {
            rSingleParticle = 0;
            for (int j = 0; j < m_nDimensions; j++) {
                rSingleParticle  += r[i][j]*r[i][j];
            }
            eKinetic += -2*alpha  + alpha*alpha * rSingleParticle;
        }
        return eKinetic;
    } else {
        return laplaceNumerical(r);
    }
}

void WaveSimple::setUseAnalyticalLaplace(bool val)
{
    useAnalytical = val;
}
