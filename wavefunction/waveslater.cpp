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
    jastrow = new Jastrow(config);

    orbitals = new Orbital*[nParticles/2];

    int shells = 0;
    if(nParticles < 3) {
        shells = 1;
    } else if(nParticles < 7) {
        shells = 2;
    } else if(nParticles < 13) {
        shells = 3;
    } else if(nParticles < 21) {
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
                if(orbital > nParticles / 2 + 1) {
                    break;
                }
            }
        }
    }
    slaterUp = new Slater(config, orbitals, true);
    slaterDown = new Slater(config, orbitals, false);
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
    for(int i = 0; i < nParticles / 2; i++) {
        orbitals[i]->setParameters(parameters);
    }
    jastrow->setParameters(parameters);
}

double WaveSlater::evaluate(vec2 r[])
{
    double normFactorial = 1;
    for(int i = 2; i < nParticles / 2; i++) {
        normFactorial *= i;
    }
    double evaluation = 1 / sqrt(normFactorial) * slaterUp->determinant(r) * slaterDown->determinant(r);
    if(m_interactionEnabled) {
        evaluation *= jastrow->evaluate(r);
    }
    return evaluation;
}

double WaveSlater::ratio(vec2 &rParticle, int particleNumber) {
    // TODO Is it necessary to calculate the inverse here?
//    slaterUp->calculateInverse();
//    slaterDown->calculateInverse();
    double theRatio = slaterUp->ratio(rParticle, particleNumber) * slaterDown->ratio(rParticle, particleNumber);
    if(m_interactionEnabled) {
        theRatio *= jastrow->ratio(rParticle, particleNumber);
    }
    return theRatio;
}

void WaveSlater::init(vec2 r[]) {
    WaveFunction::init(r);
    slaterUp->constructMatrix(r);
    std::cout << slaterUp->matrix() << std::endl;
    slaterUp->calculateInverse();
    slaterDown->constructMatrix(r);
    std::cout << slaterDown->matrix() << std::endl;
    slaterDown->calculateInverse();
    std::cout << "Doing fine" << std::endl;
}

//void WaveSlater::setPreviousMovedParticle(int particleNumber) {
//    WaveFunction::setPreviousMovedParticle(particleNumber);
//    // TODO Set last moved particle for slater or something
//}

void WaveSlater::acceptEvaluation() {
    WaveFunction::acceptEvaluation();
    slaterDown->acceptEvaluation();
    slaterUp->acceptEvaluation();
    jastrow->acceptEvaluation();
}
