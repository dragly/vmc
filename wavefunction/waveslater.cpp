#include "waveslater.h"

#include "../slater/slater.h"
#include "../jastrow/jastrow.h"
#include "../orbital/orbital.h"

#include <armadillo>

using namespace arma;

WaveSlater::WaveSlater(Config *config) :
    WaveFunction(config),
    m_interactionEnabled(config->interactionEnabled())
{
    slater = new Slater(config);
    jastrow = new Jastrow(config);

    orbitals = new Orbital*[m_nParticles/2];

    int shells = 0;
    if(m_nParticles < 3) {
        shells = 1;
    } else if(m_nParticles < 7) {
        shells = 2;
    } else if(m_nParticles < 13) {
        shells = 3;
    } else if(m_nParticles < 21) {
        shells = 4;
    } else {
        cerr << "Too many particles! Cannot calculate the number of orbitals." << endl;
        exit(918);
    }
    // generate the orbitals
    int orbital = 0;
    for(int i = 0; i < shells; i++) {
        for(int j = 0; j < shells; j++) {
            int nx = i;
            int ny = j;
            if(nx + ny < shells) {
                orbitals[orbital] = new Orbital(nx, ny, config);
                orbital++;
                if(orbital > m_nParticles / 2 + 1) {
                    break;
                }
            }
        }
    }
}

/*!
  Sets the alpha and beta parameters. This could also point to more parameters
  if needed at a later time.

  This function is specific for WaveSlater and also sets the paramters for the
  underlying orbitals

  @param parameters A pointer to an array of parameters to set.
*/
void WaveSlater::setParameters(double *parameters) {
    WaveFunction::setParameters(parameters);
    for(int i = 0; i < m_nParticles / 2; i++) {
        orbitals[i]->setParameters(parameters);
    }
    jastrow->setParameters(parameters);
}

double WaveSlater::wave(vec2 r[])
{
    if(m_interactionEnabled) {
        return slater->determinant(r, orbitals) * jastrow->evaluate(r);
    } else {
        return slater->determinant(r, orbitals);
    }
}
