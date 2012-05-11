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
    void setPreviousMovedParticle(int particleNumber);
    double ratio(vec2 &rNew, int movedParticle);
    void acceptEvaluation(int movedParticle);
    void initialize(vec2 positions[]);
    double laplace(vec2 r[], int movedParticlea);
    void gradient(const vec2 r[], int movedParticle, vec &rGradient) const;
    bool hasParticle(int particleNumber) const;
private:
    mat currentMatrix;
    mat previousMatrix;

    mat currentInverse;
    mat previousInverse;

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
