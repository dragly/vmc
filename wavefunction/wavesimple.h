#ifndef WAVESIMPLE_H
#define WAVESIMPLE_H

#include <armadillo>

#include "wavefunction.h"

using namespace std;
using namespace arma;

/*!
* \brief Used for the two particle case without interaction
  */
class WaveSimple : public WaveFunction
{
public:
    WaveSimple(Config *config);
    double evaluate(vec2 r[]);
    double laplace(vec2 r[]);
    WaveFunction* clone() {
        std::cerr << "Clone not implemented for WaveSimple" << std::endl;
        exit(978);
        return 0;
    }
private:
    bool useAnalytical;
    vec2 rPlus[];
    vec2 rMinus[];
};

#endif // WAVESIMPLE_H
