#include "waveslater.h"

#include "../slater/slater.h"
#include "../jastrow/jastrow.h"
#include "../orbital/orbital.h"

#include <armadillo>

using namespace arma;

WaveSlater::WaveSlater(Config *config) :
    WaveFunction(config)
{
    slater = new Slater(config);
    jastrow = new Jastrow(config);

    orbitals = new Orbital*[m_nParticles/2];

    // generate the orbitals
    int orbital = 0;
    int nInOrbital = 0;
    for(int i = 0; i < m_nParticles / 2; i++) {
        int nx;
        int ny;
        nx = nInOrbital;
        ny = orbital - nInOrbital;
        orbitals[i] = new Orbital(nx, ny, config);
        nInOrbital++;
        if(nInOrbital >= orbital) {
            nInOrbital = 0;
            orbital++;
        }
    }
}

void WaveSlater::setParameters(double *parameters) {
    WaveFunction::setParameters(parameters);
    for(int i = 0; i < m_nParticles / 2; i++) {
        orbitals[i]->setParameters(m_parameters);
    }
    jastrow->setParameters(m_parameters);
}

double WaveSlater::wave(const vec2 r[])
{
    // TODO implement Jastrow-factor
    return slater->determinant(r, orbitals) * jastrow->evaluate(r);
}
