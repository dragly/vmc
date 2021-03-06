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
/*!
  * \brief Used for the many-body case with interaction.
  */
class WaveSlater : public WaveFunction
{
public:
    WaveSlater(Config *config);
    double evaluate(vec2 r[]);
    void setParameters(double *parameters);
//    void setPreviousMovedParticle(int particleNumber);
    void acceptMove(int movedParticle);
    double ratio(vec2 &particlePosition, int movedParticle);
    void initialize(vec2 positions[]);
    double laplace(vec2 r[]);
    void gradient(vec2 r[], vec &rGradient);
    void rejectMove();
    WaveFunction *clone();
    Slater *slaterUp;
    Slater *slaterDown;
    Jastrow *jastrow;
    ~WaveSlater();
    void prepareGradient(vec2 &particlePosition, int movedParticle);
    void outputProperties() {
        WaveFunction::outputProperties();
        std::cout << interactionEnabled << std::endl;
    }
private:
    vec slaterUpGradient;
    vec slaterDownGradient;
    vec jastrowGradient;

    Orbital **orbitals;
    bool interactionEnabled;
};

#endif // WAVESLATER_H
