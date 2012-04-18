#include "waveslater.h"

#include "../slater/slater.h"

WaveSlater::WaveSlater(Config *config) :
    WaveFunction(config)
{
    slater = new Slater(config);
}

double WaveSlater::wave(vec2 *r)
{
    // TODO implement Jastrow-factor
    return slater->determinant(r);
}
