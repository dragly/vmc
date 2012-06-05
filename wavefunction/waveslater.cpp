#include "waveslater.h"

#include "../slater/slater.h"
#include "../jastrow/jastrow.h"
#include "../orbital/orbital.h"

#include <armadillo>

using namespace arma;

/*!
 * \brief WaveSlater::WaveSlater implements the Slater and Jastrow factor into a WaveFunction.
 * \param config A Config class with preloaded configuration
 */
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

    slaterUpGradient = zeros<vec>(nParticles * nDimensions);
    slaterDownGradient = zeros<vec>(nParticles * nDimensions);
    jastrowGradient = zeros<vec>(nParticles * nDimensions);
}

WaveSlater::~WaveSlater() {
    delete slaterUp;
    delete slaterDown;
    delete [] orbitals;
    delete jastrow;
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

/*!
 * \brief WaveSlater::evaluate Evaluates the wave function.
 * \param r The position to evaluate the wave function
 * \return The value of the wave function
 */
double WaveSlater::evaluate(vec2 r[])
{
    double normFactorial = 1;
    for(int i = 2; i < nParticles / 2; i++) {
        normFactorial *= i;
    }
//    std::cout << "Evaluation:\n " << slaterUp->matrix() << std::endl << "Norm factorial: " << normFactorial << std::endl;
    double evaluation = 1. / sqrt(normFactorial) * slaterUp->determinant(r) * slaterDown->determinant(r);
//    std::cout << "determinant:" << slaterUp->determinant(r) << std::endl << slaterUp->matrix() << std::endl;
    if(interactionEnabled) {
        evaluation *= jastrow->evaluate(r);
    }
    if(evaluation == 0) {
        std::cerr << "WARNING: Wave evaluated as zero! " << slaterUp->determinant(r) << " " << slaterDown->determinant(r) << " " << jastrow->evaluate(r) << std::endl;
    }
    return evaluation;
}

/*!
 * \brief WaveSlater::ratio Uses the Slater and Jastrow ratios to calculate the wave function ratio.
 * \param particlePosition
 * \param movedParticle
 * \return The ratio
 */
double WaveSlater::ratio(vec2 &particlePosition, int movedParticle) {
    rNew[movedParticle] = particlePosition;
//    std::cout << "after:" << std::endl;
//    for(int i = 0; i < nParticles; i++) {
//        std::cout << std::setprecision(20) << rNew[i][0] << "," << rNew[i][1] << std::endl;
//    }
    slaterUp->updateMatrix(particlePosition, movedParticle);
//    slaterUp->updateInverse(particlePosition, movedParticle);
    slaterDown->updateMatrix(particlePosition, movedParticle);
//    slaterDown->updateInverse(particlePosition, movedParticle);
    double theRatio = slaterUp->ratio(particlePosition, movedParticle) * slaterDown->ratio(particlePosition, movedParticle);
    if(interactionEnabled) {
        theRatio *= jastrow->ratio(particlePosition, movedParticle);
    }
    return theRatio;
//    return WaveFunction::ratio(rParticle, movedParticle);
}

/*!
 * \brief WaveSlater::initialize Initializes the wave function and the Slater and Jastrow factors
 * \param positions
 */
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
/*!
 * \brief WaveSlater::acceptMove Calls the acceptMove functions in Slater and Jastrow to for instance update the Slater matrix.
 * \param movedParticle
 */
void WaveSlater::acceptMove(int movedParticle) {
    WaveFunction::acceptMove(movedParticle);
    slaterUp->acceptMove(movedParticle);
    slaterDown->acceptMove(movedParticle);
    jastrow->acceptMove(movedParticle);
}

/*!
 * \brief WaveSlater::rejectMove calls the rejectMove functions in the Slater and Jastrow classes that for instance resets the Slater matrix.
 */
void WaveSlater::rejectMove() {
    WaveFunction::rejectMove();
    slaterUp->rejectMove();
    slaterDown->rejectMove();
    jastrow->rejectMove();
}

/*!
 * \brief WaveSlater::laplace calculates the laplace to wave ratio.
 * The equation used to construct this function is from eq. 16.16 in Hjorth-Jensen's Lecture Notes.
 * \param r An array of position
 * \param movedParticle The index of the particle which was last moved
 * \return Laplace to wave ratio \f$\frac{\nabla^2 \Psi}{\Psi}\f$
 */
double WaveSlater::laplace(vec2 r[]) {
    if(!useAnalyticalLaplace) {
        return laplaceNumerical(r);
    }
    double laplaceSum = 0;
    laplaceSum += slaterUp->laplace(r);
//    std::cout << "up " << slaterUp->laplace(r, movedParticle) << std::endl;
    laplaceSum += slaterDown->laplace(r);
    if(interactionEnabled) {
//    std::cout << "down " << slaterDown->laplace(r, movedParticle) << std::endl;
//        jastrow->calculateDistances(r);

        slaterUpGradient.zeros();
        slaterDownGradient.zeros();
        jastrowGradient.zeros();
        slaterUp->gradient(r, slaterUpGradient);
        slaterDown->gradient(r, slaterDownGradient);
        jastrow->gradient(r, jastrowGradient);

        laplaceSum += dot(jastrowGradient, jastrowGradient);
        laplaceSum += jastrow->laplacePartial(r);

        laplaceSum += 2 * dot(slaterUpGradient + slaterDownGradient, jastrowGradient);
    }
    return laplaceSum;
}

/*!
 * \brief WaveSlater::gradient calculates the gradient for a given configuration
 * \param r
 * \param rGradient A nParticles * nDimensions gradient where each particle's gradient is found in the i * nDimensions and (i+1) * nDimensions indices
 */
void WaveSlater::gradient(vec2 r[], vec &rGradient) {
    if(!useAnalyticalGradient) {
        return gradientNumerical(r, rGradient);
    }
    slaterUpGradient.zeros();
    slaterDownGradient.zeros();
    jastrowGradient.zeros();
    slaterUp->gradient(r, slaterUpGradient);
    slaterDown->gradient(r, slaterDownGradient);
    if(interactionEnabled) {
        jastrow->gradient(r, jastrowGradient);
    }
    rGradient = slaterUpGradient + slaterDownGradient + jastrowGradient;
//    gradientNumerical(r, movedParticle, rGradient);
}

/*!
 * \brief WaveSlater::clone is an harsh and simple function that creates a new WaveSlater object with the same config and initializes it
 * with the same positions as this object.
 * \return A WaveSlater object that will have most properties equal with this object.
 */
WaveFunction* WaveSlater::clone() {
    WaveSlater *myCopy = new WaveSlater(config);
    myCopy->setUseAnalyticalGradient(this->useAnalyticalGradient);
    myCopy->setUseAnalyticalLaplace(this->useAnalyticalLaplace);
    myCopy->interactionEnabled = this->interactionEnabled;
    myCopy->setParameters(parameters);
    myCopy->previousEvaluation = this->previousEvaluation;
    myCopy->currentEvaluation = this->currentEvaluation;
    myCopy->nParticles = this->nParticles;
    myCopy->nDimensions = this->nDimensions;
    for(int i = 0; i < nParticles; i++) {
        myCopy->rPlus[i] = this->rPlus[i];
        myCopy->rMinus[i] = this->rMinus[i];
        myCopy->rNew[i] = this->rNew[i];
        myCopy->rOld[i] = this->rOld[i];
    }
    myCopy->slaterUpGradient = this->slaterUpGradient;
    myCopy->slaterDownGradient = this->slaterDownGradient;
    myCopy->jastrowGradient = this->jastrowGradient;
//    myCopy->initialize(rOld);
    return myCopy;
}

/*!
 * \brief WaveSlater::prepareGradient updates the Slater matrix and inverse to prepare for gradient calculations
 * \param particlePosition The new position for the moved particle
 * \param movedParticle The index of the moved particle
 */
void WaveSlater::prepareGradient(vec2 &particlePosition, int movedParticle) {
    rNew[movedParticle] = particlePosition;
    slaterUp->updateMatrix(particlePosition, movedParticle);
    slaterUp->updateInverse(particlePosition, movedParticle);
//    slaterUp->calculateInverseNumerically();
    slaterDown->updateMatrix(particlePosition, movedParticle);
    slaterDown->updateInverse(particlePosition, movedParticle);
//    slaterDown->calculateInverseNumerically();
//    jastrow->calculateDistances(rNew);
}
