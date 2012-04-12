#ifndef SLATER_H
#define SLATER_H

#include <armadillo>
using namespace arma;

class Config;

class Slater
{
public:
    Slater(Config *config);

    ~Slater();
private:
    mat *matrix;

    int m_nDimensions;
    int m_nParticles;
};

#endif // SLATER_H
