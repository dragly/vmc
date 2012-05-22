#include "slater.h"

#include "../config.h"
#include "../orbital/orbital.h"

#include <armadillo>
#include <iomanip>
#include <iostream>

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
    if(spinUp) {
        particleIndexOffset = 0;
    } else {
        particleIndexOffset = nParticles / 2;
    }
}

void Slater::initialize(vec2 positions[])
{
    constructMatrix(positions);
    calculateInverseNumerically();
    currentRatio = 1;
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
    currentMatrix = previousMatrix;
    // TODO Remove this
//    previousMatrix = randu<mat>(nParticles / 2, nParticles / 2);
}

double Slater::determinant(vec2 r[]) {
    constructMatrix(r);
    double determ = det(previousMatrix);
    return determ;
}

/*!

  Note: You need to run ratio() before calling this function.
  */
void Slater::calculateInverse(int movedParticle)
{
    if(hasParticle(movedParticle)) {
        int localParticle = movedParticle - particleIndexOffset;
        int p = localParticle;
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
//    calculateInverseNumerically(); // TODO revert this to analytical solution of inverse after fixing it
}

int nExceptions = 0;
void Slater::calculateInverseNumerically() {
    try {
        previousInverse = inv(previousMatrix);
        currentInverse = previousInverse;
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

void Slater::setPreviousMovedParticle(int movedParticle)
{
    previousMovedParticle = movedParticle;
}

void Slater::updateMatrix(vec2 &particlePosition, int movedParticle) {
    if(hasParticle(movedParticle)) {
        int localParticle = movedParticle - particleIndexOffset;
        for(int j = 0; j < nParticles / 2; j++) {
            currentMatrix.at(localParticle,j) = orbitals[j]->evaluate(particlePosition);
        }
    }
}

/*!
  Note: The first half of the particles have spin up, while the others have spin down.
  */
double Slater::ratio(vec2 &particlePosition, int movedParticle)
{
    if(hasParticle(movedParticle)) {
        int localParticle = movedParticle - particleIndexOffset;
        updateMatrix(particlePosition, movedParticle);
        movedParticle = movedParticle - particleIndexOffset;
        double R = 0;
        for(int i = 0; i < nParticles / 2; i++) {
            R += currentMatrix.at(localParticle,i) * previousInverse.at(i, localParticle);
        }
        currentRatio = R;
        return R;
    } else {
        currentRatio = 1;
        return 1;
    }
}

bool Slater::hasParticle(int particleNumber) const {
    return (spinUp && particleNumber < nParticles / 2) || (!spinUp && particleNumber >= nParticles / 2);
}

void Slater::acceptMove(int movedParticle)
{
    previousMatrix = currentMatrix;
    previousRatio = currentRatio;
    calculateInverse(movedParticle);
    previousInverse = currentInverse;
}

void Slater::rejectMove() {
    currentMatrix = previousMatrix;
    // TODO might not need to copy back the inverse
    currentInverse = previousInverse;
    currentRatio = previousRatio;
}

void Slater::gradient(vec2 r[], vec &rGradient) {
    rGradient.zeros();
    for(int a = 0; a < nParticles; a++) {
        // TODO we are now recalculating the gradient for all particles, this could be avoided
        int movedParticle = a;
        if(hasParticle(movedParticle)) {
            int localParticle = movedParticle - particleIndexOffset;
            for(int j = 0; j < nParticles / 2; j++) {
                orbitals[j]->gradient(r[movedParticle], orbitalGradient);
                double inverseValue = previousInverse(j,localParticle);
                for(int l = 0; l < nDimensions; l++) {
                    int rGradientIndex = movedParticle * nDimensions + l;
                    rGradient[rGradientIndex] += orbitalGradient[l] * inverseValue;
                }
            }
        }
    }
}

double Slater::laplace(vec2 r[])
{
    // TODO we are now recalculating the laplace for all particles, this could be avoided
    double laplaceResult = 0;
    for(int a = 0; a < nParticles; a++) {
        // TODO we are now recalculating the gradient for all particles, this could be avoided
        int movedParticle = a;
        if(hasParticle(movedParticle)) {
            int localParticle = movedParticle - particleIndexOffset;
//            std::cout << "a" << localParticle << std::endl;
            for(int j = 0; j < nParticles / 2; j++) {
                laplaceResult += orbitals[j]->laplace(r[movedParticle]) * currentInverse(j,localParticle);
//                std::cout << orbitals[j]->laplace(r[movedParticle]) << std::endl;
//                std::cout << "inverse: " << previousInverse(j,localParticle) << std::endl;
            }
        }
    }
    return laplaceResult;
}

mat Slater::inverse() {
    return currentInverse;
}

mat Slater::matrix() {
    return currentMatrix;
}

Slater::~Slater() {
}
