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
    void calculateInverse(int movedParticle);
    void calculateInverseNumerically();
    void setPreviousMovedParticle(int movedParticle);
    double ratio(vec2 &particlePosition, int movedParticle);
    void acceptEvaluation(int movedParticle);
    void initialize(vec2 positions[]);
    double laplace(vec2 r[], int movedParticlea);
    void gradient(vec2 r[], int movedParticle, vec &rGradient);
    bool hasParticle(int particleNumber) const;
    void refuseEvaluation();
    void updateMatrix(vec2 &particlePosition, int movedParticle);
private:
    mat currentMatrix;
    mat previousMatrix;

    mat currentInverse;
    mat previousInverse;

    vec2 orbitalGradient;

    int nDimensions;
    int nParticles;

    int previousMovedParticle;

    double previousRatio;
    double currentRatio;

    int particleIndexOffset;

    Orbital **orbitals;

    bool spinUp;
};

#endif // SLATER_H
