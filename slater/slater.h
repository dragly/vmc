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
    void initialize(vec2 r[]);
    double laplace(const vec2 &r);
    void gradient(const vec2 r[], vec2 &rGradient) const;
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

    Orbital **orbitals;

    bool spinUp;
};

#endif // SLATER_H
