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
    previousMatrix = zeros<mat>(nParticles / 2, nParticles / 2);
    currentMatrix = zeros<mat>(nParticles / 2, nParticles / 2);
    previousInverse = zeros<mat>(nParticles / 2, nParticles / 2);
    currentInverse = zeros<mat>(nParticles / 2, nParticles / 2);
}

void Slater::initialize(vec2 r[])
{
    constructMatrix(r);
    calculateInverseNumerically();
    currentRatio = 1;
}

double Slater::laplace(const vec2 &r)
{
    return 0;
}

/*!
  Note: The first half of the particles have spin up, while the others have spin down.
  */
void Slater::constructMatrix(vec2 r[]) {
    for(int i = 0; i < nParticles / 2; i++) {
        for(int j = 0; j < nParticles / 2; j++) {
            int index = -1;
            if(spinUp) {
                index = i;
            } else {
                index = i + nParticles / 2;
            }
            previousMatrix(i,j) = orbitals[j]->evaluate(r[index]);
        }
    }
}

double Slater::determinant(vec2 r[]) {
    constructMatrix(r);
    // the slater determinant can be divided into two parts multiplied together
    double determ = det(previousMatrix);
    //    double detUp = 1;
    //    double detDown = 2;
    return determ;
}


void Slater::calculateInverse(int movedParticle)
{
    int p = movedParticle;
    vec S = zeros<vec>(nParticles / 2);
    for(int j = 0; j < nParticles / 2; j++) {
        for(int l = 0; l < nParticles / 2; l++) {
            S[j] += currentMatrix.at(p, l) * previousInverse.at(l,j);
        }
    }
    // TODO split into two for loops (if efficient)
    for(int i = 0; i < nParticles / 2; i++) {
        for(int j = 0; j < nParticles / 2; j++) {
            if(p != j) {
                currentInverse.at(i,j) = previousInverse.at(i,j) - previousInverse.at(i,p) * S.at(j) / currentRatio;
            } else {
                currentInverse.at(i,j) = previousInverse.at(i,j) / currentRatio;
            }
        }
    }
}

int nExceptions = 0;
void Slater::calculateInverseNumerically() {
    try {
        previousInverse = inv(previousMatrix);
    } catch(std::runtime_error &e) {
        std::cout << matrix() << std::endl;
        std::cout << "Caught inverse exception! Setting inverse to large random matrix." << std::endl;
        previousInverse = randu<mat>(nParticles / 2, nParticles / 2) * 99999999;
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
    bool hasParticle = (spinUp && movedParticle < nParticles / 2) || (!spinUp && movedParticle >= nParticles / 2);
    if(hasParticle) {
        for(int i = 0; i < nParticles / 2; i++) {
            currentMatrix.at(movedParticle,i) = orbitals[i]->evaluate(rNew);
        }
        if(spinUp) {
            movedParticle = movedParticle;
        } else {
            movedParticle = movedParticle - nParticles / 2;
        }
        double R = 0;
        for(int i = 0; i < nParticles / 2; i++) {
            R += currentMatrix.at(movedParticle,i) * previousInverse.at(i, movedParticle);
        }
        currentRatio = R;
        return R;
    } else {
        currentRatio = 1;
        return 1;
    }
}

void Slater::acceptEvaluation(int movedParticle)
{
    calculateInverse(movedParticle);
    previousMatrix = currentMatrix;
    previousInverse = currentInverse;
    previousRatio = currentRatio;
}

void Slater::gradient(const vec2 r[], vec2 &rGradient) const {

}

mat Slater::inverse() {
    return previousInverse;
}

mat Slater::matrix() {
    return previousMatrix;
}

Slater::~Slater() {
}
