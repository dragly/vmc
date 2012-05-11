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
    std::cout << "Creating with " << nParticles << std::endl;
    jastrow = new Jastrow(config);

    orbitals = new Orbital*[nParticles/2];

    int shells = 0;
    if(nParticles <= 2) {
        shells = 1;
    } else if(nParticles <= 6) {
        shells = 2;
    } else if(nParticles <= 12) {
        shells = 3;
    } else if(nParticles <= 20) {
        shells = 4;
    } else {
        cerr << "Too many particles! Cannot calculate the number of orbitals for more than 20 particles." << endl;
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
//    std::cout << "Evaluation:\n " << slaterUp->matrix() << std::endl << "Norm factorial: " << normFactorial << std::endl;
    double evaluation = 1 / sqrt(normFactorial) * slaterUp->determinant(r) * slaterDown->determinant(r);
//    std::cout << "determinant:" << slaterUp->determinant(r) << std::endl << slaterUp->matrix() << std::endl;
    if(m_interactionEnabled) {
        evaluation *= jastrow->evaluate(r);
    }
    return evaluation;
}

double WaveSlater::ratio(vec2 &particlePosition, int particleNumber) {
    double theRatio = slaterUp->ratio(particlePosition, particleNumber) * slaterDown->ratio(particlePosition, particleNumber);
    if(m_interactionEnabled) {
        theRatio *= jastrow->ratio(particlePosition, particleNumber);
    }
    return theRatio;
//    return WaveFunction::ratio(rParticle, particleNumber);
}

void WaveSlater::initialize(vec2 positions[]) {
    slaterUp->initialize(positions);
    slaterDown->initialize(positions);
    WaveFunction::initialize(positions);
    std::cout << "Slater up init:\n" << slaterUp->matrix() << std::endl << slaterUp->inverse() << std::endl;
}

//void WaveSlater::setPreviousMovedParticle(int particleNumber) {
//    WaveFunction::setPreviousMovedParticle(particleNumber);
//    // TODO Set last moved particle for slater or something
//}

void WaveSlater::acceptEvaluation(int movedParticle) {
    WaveFunction::acceptEvaluation(movedParticle);
    slaterDown->acceptEvaluation(movedParticle);
    slaterUp->acceptEvaluation(movedParticle);
    jastrow->acceptEvaluation(movedParticle);
}

double WaveSlater::laplace(vec2 r[], int movedParticle) {
    double laplaceSum = 0;
    laplaceSum += slaterUp->laplace(r, movedParticle);
    laplaceSum += slaterDown->laplace(r, movedParticle);
//    laplaceSum += jastrow->laplace(r, movedParticle);
    return laplaceSum;
}

void WaveSlater::gradient(vec2 r[], int particleNumber, vec &rGradient) {
    std::cout << "rGrad size " << rGradient.size() << " " << nParticles * nDimensions << std::endl;
    vec slaterUpGradient = zeros<vec>(nParticles * nDimensions);
    vec slaterDownGradient = zeros<vec>(nParticles * nDimensions);
    vec jastrowGradient = zeros<vec>(nParticles * nDimensions);
    slaterUp->gradient(r, particleNumber, slaterUpGradient);
    slaterDown->gradient(r, particleNumber, slaterDownGradient);
    std::cout << "Slater Up\n" << slaterUpGradient << std::endl;
    std::cout << "Slater Down\n" << slaterDownGradient << std::endl;
//    jastrow->gradient(r, jastrowGradient);
    rGradient = slaterUpGradient + slaterDownGradient + jastrowGradient;
    std::cout << "Total gradient\n" << rGradient << std::endl;
//    gradientNumerical(r, particleNumber, rGradient);
}
