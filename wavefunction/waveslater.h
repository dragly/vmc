#ifndef WAVESLATER_H
#define WAVESLATER_H

#include "wavefunction.h"

#include <armadillo>
using namespace std;
using namespace arma;

#include "../config.h"

class Slater;
class Jastrow;
class Orbital;

/*!
  Wavefunction class that uses the Slater determinant to handle a general
  number of particles.
*/
class WaveSlater : public WaveFunction
{
public:
    WaveSlater(Config *config);
    double evaluate(vec2 r[]);
    void setParameters(double *parameters);
//    void setPreviousMovedParticle(int particleNumber);
    void acceptEvaluation(int movedParticle);
    double ratio(vec2 &particlePosition, int movedParticle);
    void initialize(vec2 positions[]);
    double laplace(vec2 r[], int movedParticle);
    void gradient(vec2 r[], int movedParticle, vec &rGradient);
    void refuseEvaluation();
private:
    Slater *slaterDown;
    Slater *slaterUp;
    Jastrow *jastrow;

    Orbital **orbitals;
    bool interactionEnabled;
};

#endif // WAVESLATER_H
