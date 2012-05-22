#include "waveideal.h"
#include <math.h>

WaveIdeal::WaveIdeal(Config *config) :
    WaveFunction(config)
{
}

double WaveIdeal::evaluate(vec2 r[])
{
    double alpha = parameters[0];
    double beta = parameters[1];
    int i, j;
    double wf, argument, r_single_particle;

    argument = wf = 0;
    for (i = 0; i < nParticles; i++) {
        r_single_particle = 0;
        for (j = 0; j < nDimensions; j++) {
            r_single_particle  += r[i][j]*r[i][j];
        }
        argument += r_single_particle;
    }
    double vec[2];
    vec[0] = r[0][0]-r[1][0];
    vec[1] = r[0][1]-r[1][1];
    double r12 = sqrt(vec[0]*vec[0]+vec[1]*vec[1]);
    double a = 1;
    wf = exp(-(argument*alpha) / 2);
    double jastrowArgument = (a * r12) / (1 + beta * r12);
    wf *= exp(jastrowArgument);
    return wf;
}

double WaveIdeal::laplace(vec2 r[])
{
    double alpha = parameters[0];
    double beta = parameters[1];
    if(useAnalyticalLaplace) {
        double omega = 1;
        double aconst = 1;
        double r12vec[2];
        r12vec[0] = r[0][0]-r[1][0];
        r12vec[1] = r[0][1]-r[1][1];
        double r12 = sqrt(r12vec[0]*r12vec[0]+r12vec[1]*r12vec[1]);
        double dotProduct = r[0][0]*r[1][0]+r[0][1]*r[1][1];
        double denom = (1 + beta*r12);

        double crossProd = 0;
        double laplaceE = 0;
        for (int i = 0; i < nParticles; i++) {
            double rSquared = 0;
            for (int j = 0; j < nDimensions; j++) {
                rSquared  += r[i][j]*r[i][j];
            }
            laplaceE += - 2*omega*alpha + omega*omega*alpha*alpha * rSquared;
            crossProd += - omega * alpha * (rSquared - dotProduct) * aconst / (r12 * denom*denom);
        }
        double laplaceC = 2 *(aconst/(denom*denom)) * (aconst/(denom*denom) - 2*beta/denom + 1/r12);
        double eKinetic = laplaceE + laplaceC + 2*crossProd;
        return eKinetic;
    } else {
        return laplaceNumerical(r);
    }
}
