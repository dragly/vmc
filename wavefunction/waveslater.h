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
    void acceptEvaluation();
    double ratio(vec2 &rParticle, int particleNumber);
    void init(vec2 r[]);
private:
    Slater *slaterDown;
    Slater *slaterUp;
    Jastrow *jastrow;

    Orbital **orbitals;
    bool m_interactionEnabled;
};

#endif // WAVESLATER_H
