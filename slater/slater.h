#ifndef SLATER_H
#define SLATER_H

#include <armadillo>
using namespace arma;

class Config;
class Orbital;

class Slater
{
public:
    Slater(Config *config);

    ~Slater();
    double determinant(vec2 *r);
    void constructMatrix(vec2 *r);
private:
    mat matrixUp;
    mat matrixDown;

    int nDimensions;
    int nParticles;

    Orbital **orbitals;
};

#endif // SLATER_H
