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

    mat matrix();
    mat inverse();
    void calculateInverse();
    void setPreviousMovedParticle(int particleNumber);
    double ratio(vec2 &rNew, int movedParticle);
    void acceptEvaluation();
private:
    mat matrixNew;
    mat matrixOld;

    mat inverseNew;
    mat inverseOld;

    int nDimensions;
    int nParticles;

    int previousMovedParticle;

    Orbital **orbitals;

    bool spinUp;
};

#endif // SLATER_H
