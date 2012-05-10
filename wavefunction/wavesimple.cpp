#include "wavesimple.h"
#include "../matrix.h"
#include "../utils.h"
#include <math.h>
#include <iostream>
using namespace std;

WaveSimple::WaveSimple(Config *config)  :
    WaveFunction(config)
{
}

double WaveSimple::evaluate(vec2 r[])
{
    double alpha = parameters[0];
    int i, j;
    double wf, argument, rSingleParticle;

    argument = wf = 0;
    for (i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for (j = 0; j < nDimensions; j++) {
            rSingleParticle  += r[i][j]*r[i][j];
        }
        argument += rSingleParticle;
    }
    wf = exp(-(argument*alpha) / 2) ;
    return wf;
}

double WaveSimple::laplace(vec2 r[])
{
    double alpha = parameters[0];
    if(useAnalyticalLaplace) {
        double eKinetic = 0;
        double rSingleParticle;
        for (int i = 0; i < nParticles; i++) {
            rSingleParticle = 0;
            for (int j = 0; j < nDimensions; j++) {
                rSingleParticle  += r[i][j]*r[i][j];
            }
            eKinetic += -2*alpha  + alpha*alpha * rSingleParticle;
        }
        return eKinetic;
    } else {
        return laplaceNumerical(r);
    }
}
