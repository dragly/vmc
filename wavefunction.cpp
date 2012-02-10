#include "wavefunction.h"

#include "matrix.h"
#include "utils.h"

WaveFunction::WaveFunction(int nParticles, int dimensions) :
    nParticles(nParticles),
    dimensions(dimensions)
{
    // allocate matrices which contain the position of the particles
    // the function matrix is defined in the progam library
    rPlus = (double **) matrix( nParticles, dimensions, sizeof(double));
    rMinus = (double **) matrix( nParticles, dimensions, sizeof(double));
}

double WaveFunction::laplaceNumerical(double **r)
{
    double eKinetic = 0;
    double wfold = wave(r);
    for (int i = 0; i < nParticles; i++) {
        for (int j=0; j < dimensions; j++) {
            rPlus[i][j] = rMinus[i][j] = r[i][j];
        }
    }
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < dimensions; j++) {
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

void WaveFunction::setParameters(double alpha, double beta)
{
    this->alpha = alpha;
    this->beta = beta;
}


WaveFunction::~WaveFunction()
{
    free_matrix((void **) rPlus); // free memory
    free_matrix((void **) rMinus);
}
