#ifndef JASTROW_H
#define JASTROW_H

#include <armadillo>

using namespace arma;

class Config;

class Jastrow
{
public:
    Jastrow(Config * config);

    double evaluate(vec2 r[]);
    double ratio(vec2 &rNew, int movedParticle);
    void setParameters(double *parameters);
    void calculateDistances(vec2 r[]);
    void initialize(vec2 positions[]);
    void acceptEvaluation(int movedParticle);
    double argument(int i, int j, mat &distances);

    void gradient(vec2 r[], int movedParticle, vec &rGradient);
    void refuseEvaluation();
    ~Jastrow();
    double laplace(vec2 r[], int movedParticle);
private:
    int nParticles;
    int nDimensions;
    double beta;
    mat a;

    mat distancesOld;
    mat distancesNew;
    mat jastrowArgumentsOld;
    mat jastrowArgumentsNew;

    vec2 *rOld;
    vec2 *rNew;
};

#endif // JASTROW_H
