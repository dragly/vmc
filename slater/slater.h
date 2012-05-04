#ifndef SLATER_H
#define SLATER_H

#include <armadillo>
using namespace arma;

class Config;
class Orbital;

class Slater
{
public:
    Slater(Config *config, Orbital *orbitals[], bool spinUp_);

    ~Slater();
    double determinant(vec2 r[]);
    void constructMatrix(vec2 r[]);
private:
    mat matrixUpNew;
    mat matrixDownNew;
    mat matrixOld;
    mat matrixDownOld;

    int nDimensions;
    int nParticles;

    Orbital **orbitals;

    bool spinUp;
};

#endif // SLATER_H
