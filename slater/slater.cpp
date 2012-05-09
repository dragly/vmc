#include "slater.h"

#include "../config.h"
#include "../orbital/orbital.h"

#include <armadillo>

using namespace arma;

Slater::Slater(Config *config, Orbital* orbitals[], bool spinUp_) :
    nDimensions(config->nDimensions()),
    nParticles(config->nParticles()),
    orbitals(orbitals),
    spinUp(spinUp_)
{
    // the slater determinant can de separated into two parts, one for spin up and the other for spin down
    matrixOld = zeros<mat>(nParticles / 2, nParticles / 2);
    matrixNew = zeros<mat>(nParticles / 2, nParticles / 2);
    //    matrixDownOld = zeros<mat>(nParticles / 2, nParticles / 2);
}

/*!
  Note: The first half of the particles have spin up, while the others have spin down.
  */
void Slater::constructMatrix(vec2 r[]) {
    for(int i = 0; i < nParticles / 2; i++) {
        for(int j = 0; j < nParticles / 2; j++) {
            if(spinUp) {
                matrixOld(i,j) = orbitals[j]->evaluate(r[i]);
            } else {
                matrixOld(i,j) = orbitals[j]->evaluate(r[i + nParticles / 2]);
            }
        }
    }
}

double Slater::determinant(vec2 r[]) {
    constructMatrix(r);
    // the slater determinant can be divided into two parts multiplied together
    double determ = det(matrixOld);
    //    double detUp = 1;
    //    double detDown = 2;
    return determ;
}

int nExceptions = 0;

void Slater::calculateInverse() {
    try {
        inverseOld = inv(matrixOld);
    } catch(std::runtime_error &e) {
        std::cout << matrix() << std::endl;
        std::cout << "Caught inverse exception! Setting inverse to large random matrix." << std::endl;
        inverseOld = randu<mat>(nParticles / 2, nParticles / 2) * 99999999;
        std::cout << inverseOld << std::endl;
        nExceptions++;
        if(nExceptions > 10) {
            std::cout << "Too may exceptions. Quitting!" << std::endl;
            exit(1337);
        }
    }
}

void Slater::setPreviousMovedParticle(int particleNumber)
{
    previousMovedParticle = particleNumber;
}

/*!
  Note: The first half of the particles have spin up, while the others have spin down.
  */
double Slater::ratio(vec2 &rNew, int movedParticle)
{
    bool hasParticle = (spinUp && movedParticle < nParticles / 2) || (!spinUp && movedParticle > nParticles / 2);
    if(hasParticle) {
        if(spinUp) {
            movedParticle = movedParticle;
        } else {
            movedParticle = movedParticle - nParticles / 2;
        }
        double R = 0;
        for(int i = 0; i < nParticles / 2; i++) {
            matrixNew.at(movedParticle,i) = orbitals[i]->evaluate(rNew);
            R += matrixNew.at(movedParticle,i) * inverseOld.at(i, movedParticle);
        }
        return R;
    } else {
        return 1;
    }
}

void Slater::acceptEvaluation()
{
    matrixOld = matrixNew;
}

mat Slater::inverse() {
    return inverseOld;
}

mat Slater::matrix() {
    return matrixOld;
}

Slater::~Slater() {
}
