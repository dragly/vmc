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
    double determinant(vec2 r[], Orbital *orbitals[]);
    void constructMatrix(vec2 r[], Orbital *orbitals[]);
private:
    mat matrixUp;
    mat matrixDown;

<<<<<<< HEAD
    int nDimensions;
    int nParticles;

    Orbital **orbitals;
=======
    int m_nDimensions;
    int m_nParticles;
>>>>>>> b4efb06c7f2fbc3cdd14718f48638ca10eb4022d
};

#endif // SLATER_H
