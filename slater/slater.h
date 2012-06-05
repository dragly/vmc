#ifndef SLATER_H
#define SLATER_H

#include <armadillo>
using namespace arma;

class Config;
class Orbital;

/*!
  * \brief Handles the single particle wave functions in a matrix
  *
  * Also calculates the inverse and the determinant of this matrix.
  */
class Slater
{
public:
    Slater(Config *config, Orbital *orbitals[], bool spinUp_);

    ~Slater();
    double determinant(vec2 r[]);
    void constructMatrix(vec2 r[]);

    mat matrix();
    mat inverse();
    void updateInverse(vec2 &particlePosition, int movedParticle);
    void calculateInverseNumerically();
    void setPreviousMovedParticle(int movedParticle);
    double ratio(vec2 &particlePosition, int movedParticle);
    void acceptMove(int movedParticle);
    void initialize(vec2 positions[]);
    double laplace(vec2 r[]);
    void gradient(vec2 r[], vec &rGradient);
    bool hasParticle(int particleNumber) const;
    void rejectMove();
    void updateMatrix(vec2 &particlePosition, int movedParticle);
private:
    vec2 *rOld;
    vec2 *rNew;

    mat currentMatrix;
    mat previousMatrix;

    mat currentInverse;
    mat previousInverse;

    vec2 orbitalGradient;

    int nDimensions;
    int nParticles;

    int previousMovedParticle;

    int particleIndexOffset;

    Orbital **orbitals;

    bool spinUp;
};

#endif // SLATER_H
