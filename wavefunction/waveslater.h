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
    double ratio(vec2 &particlePosition, int particleNumber);
    void initialize(vec2 positions[]);
    double laplace(vec2 r[]);
    void gradient(vec2 r[], vec2 &rGradient);
private:
    Slater *slaterDown;
    Slater *slaterUp;
    Jastrow *jastrow;

    Orbital **orbitals;
    bool m_interactionEnabled;
};

#endif // WAVESLATER_H
