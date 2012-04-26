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
    void setParameters(double *parameters);
private:
    int m_nParticles;
    double m_beta;
    mat a;
};

#endif // JASTROW_H
