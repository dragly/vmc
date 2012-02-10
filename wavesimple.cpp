#include "wavesimple.h"
#include "matrix.h"
#include "utils.h"
#include <math.h>
#include <iostream>
using namespace std;

WaveSimple::WaveSimple(int number_particles, int dimension)  :
    WaveFunction(),
    number_particles(number_particles),
    dimension(dimension),
    useAnalytical(false),
    nExp(100000)
{
    // allocate matrices which contain the position of the particles
    // the function matrix is defined in the progam library
    r_plus = (double **) matrix( number_particles, dimension, sizeof(double));
    r_minus = (double **) matrix( number_particles, dimension, sizeof(double));
    preExp = new double[nExp];
    aExp = -6;
    bExp = -3;
    ratioExp = nExp / (bExp - aExp);
    ratioExpInverse = (bExp - aExp) / nExp;
    nExpInverse = 1. / nExp;

    for(int i = 0; i < nExp; i++) {
        double x = aExp + i * ratioExpInverse;
        preExp[i] = exp(x);
    }
}

double WaveSimple::myExp(double x) {
    int lowIndex = (int)((x - aExp) * ratioExp);
    if(lowIndex < 0 || lowIndex > nExp) {
//        cerr << "Index error. x value was " << x << endl;
        return exp(x);
    } else {
        return preExp[lowIndex];
    }
}

double WaveSimple::wave(double **r)
{
    int i, j;
    double wf, argument, rSingleParticle;

    argument = wf = 0;
    for (i = 0; i < number_particles; i++) {
        rSingleParticle = 0;
        for (j = 0; j < dimension; j++) {
            rSingleParticle  += r[i][j]*r[i][j];
        }
        argument += rSingleParticle;
    }
    wf = myExp(-(argument*alpha) / 2) ;
    return wf;
}

double WaveSimple::laplace(double **r)
{
    double eKinetic = 0;
    if(useAnalytical) {
        int i, j;
        double argument;
        double rSingleParticle;

        argument = 0;
        double prefix = 0;
        for (i = 0; i < number_particles; i++) {
            rSingleParticle = 0;
            for (j = 0; j < dimension; j++) {
                rSingleParticle  += r[i][j]*r[i][j];
            }
            argument += rSingleParticle;
            prefix += rSingleParticle;
        }
        eKinetic = prefix * alpha * exp(-(argument*alpha) / 2) ;
    } else {
        for (int i = 0; i < number_particles; i++) {
            for (int j=0; j < dimension; j++) {
                r_plus[i][j] = r_minus[i][j] = r[i][j];
            }
        }
        for (int i = 0; i < number_particles; i++) {
            for (int j = 0; j < dimension; j++) {
                r_plus[i][j] = r[i][j]+h;
                r_minus[i][j] = r[i][j]-h;
                double wfminus = wave(r_minus);
                double wfplus  = wave(r_plus);
                double wfold = wave(r);
                eKinetic += (wfminus+wfplus-2*wfold);
                r_plus[i][j] = r[i][j];
                r_minus[i][j] = r[i][j];
            }
        }
    }
    return eKinetic;
}

void WaveSimple::setUseAnalytical(bool val)
{
    useAnalytical = val;
}

WaveSimple::~WaveSimple()
{
    free_matrix((void **) r_plus); // free memory
    free_matrix((void **) r_minus);
}
