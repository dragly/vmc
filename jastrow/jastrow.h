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
    double ratio(vec2 &rNew, int particleNumber);
    void setParameters(double *parameters);
    void calculateDistances(vec2 r[]);
    void acceptEvaluation();
    double argument(int i, int j, mat &distances);
    ~Jastrow();
private:
    int nParticles;
    double m_beta;
    mat a;

    mat distancesOld;
    mat distancesNew;
    mat jastrowArgumentsOld;
    mat jastrowArgumentsNew;

    vec2 *rOld;
    vec2 *rNew;
};

#endif // JASTROW_H
