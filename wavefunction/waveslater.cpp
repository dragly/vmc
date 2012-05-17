#include "waveslater.h"

#include "../slater/slater.h"
#include "../jastrow/jastrow.h"
#include "../orbital/orbital.h"

#include <armadillo>

using namespace arma;

WaveSlater::WaveSlater(Config *config) :
    WaveFunction(config),
    interactionEnabled(config->interactionEnabled())
{
//    std::cout << "Creating with " << nParticles << std::endl;
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
    if(interactionEnabled) {
        evaluation *= jastrow->evaluate(r);
    }
    return evaluation;
}

double WaveSlater::ratio(vec2 &particlePosition, int movedParticle) {
    rNew[movedParticle] = particlePosition;
//    std::cout << "after:" << std::endl;
//    for(int i = 0; i < nParticles; i++) {
//        std::cout << std::setprecision(20) << rNew[i][0] << "," << rNew[i][1] << std::endl;
//    }
    double theRatio = slaterUp->ratio(particlePosition, movedParticle) * slaterDown->ratio(particlePosition, movedParticle);
    if(interactionEnabled) {
        theRatio *= jastrow->ratio(particlePosition, movedParticle);
    }
    return theRatio;
//    return WaveFunction::ratio(rParticle, movedParticle);
}

void WaveSlater::initialize(vec2 positions[]) {
    slaterUp->initialize(positions);
    slaterDown->initialize(positions);
    jastrow->initialize(positions);
    WaveFunction::initialize(positions);
}

//void WaveSlater::setPreviousMovedParticle(int movedParticle) {
//    WaveFunction::setPreviousMovedParticle(movedParticle);
//    // TODO Set last moved particle for slater or something
//}
// TODO update matrices with new values for moved particle
void WaveSlater::acceptEvaluation(int movedParticle) {
    WaveFunction::acceptEvaluation(movedParticle);
    slaterUp->acceptEvaluation(movedParticle);
    slaterDown->acceptEvaluation(movedParticle);
    jastrow->acceptEvaluation(movedParticle);
}

void WaveSlater::refuseEvaluation() {
    WaveFunction::refuseEvaluation();
    slaterUp->refuseEvaluation();
    slaterDown->refuseEvaluation();
    jastrow->refuseEvaluation();
}

/*!
 * \brief WaveSlater::laplace calculates the laplace to wave ratio.
 * The equation used to construct this function is from eq. 16.16 in Hjorth-Jensen's Lecture Notes.
 * \param r An array of position
 * \param movedParticle The index of the particle which was last moved
 * \return Laplace to wave ratio \f$\frac{\del^2 \Psi}{\Psi}\f$
 */
double WaveSlater::laplace(vec2 r[], int movedParticle) {
    if(!useAnalyticalLaplace) {
        return laplaceNumerical(r);
    }
    double laplaceSum = 0;
    laplaceSum += slaterUp->laplace(r, movedParticle);
//    std::cout << "up " << slaterUp->laplace(r, movedParticle) << std::endl;
    laplaceSum += slaterDown->laplace(r, movedParticle);
    if(interactionEnabled) {
//    std::cout << "down " << slaterDown->laplace(r, movedParticle) << std::endl;
//        jastrow->calculateDistances(r);

        vec slaterUpGradient = zeros<vec>(nParticles * nDimensions);
        vec slaterDownGradient = zeros<vec>(nParticles * nDimensions);
        vec jastrowGradient = zeros<vec>(nParticles * nDimensions);
        slaterUp->gradient(r, movedParticle, slaterUpGradient);
        slaterDown->gradient(r, movedParticle, slaterDownGradient);
        jastrow->gradient(r, movedParticle, jastrowGradient);

        laplaceSum += dot(jastrowGradient, jastrowGradient);
        laplaceSum += jastrow->laplacePartial(r, movedParticle);

        laplaceSum += 2 * dot(slaterUpGradient + slaterDownGradient, jastrowGradient);
    }
    return laplaceSum;
}

void WaveSlater::gradient(vec2 r[], int movedParticle, vec &rGradient) {
    if(!useAnalyticalGradient) {
        return gradientNumerical(r, rGradient);
    }
    vec slaterUpGradient = zeros<vec>(nParticles * nDimensions);
    vec slaterDownGradient = zeros<vec>(nParticles * nDimensions);
    vec jastrowGradient = zeros<vec>(nParticles * nDimensions);
    slaterUp->gradient(r, movedParticle, slaterUpGradient);
    slaterDown->gradient(r, movedParticle, slaterDownGradient);
    if(interactionEnabled) {
        jastrow->gradient(r, movedParticle, jastrowGradient);
    }
    rGradient = slaterUpGradient + slaterDownGradient + jastrowGradient;
//    gradientNumerical(r, movedParticle, rGradient);
}
