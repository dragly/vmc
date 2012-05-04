#include "slater.h"

#include "../config.h"
#include "../orbital/orbital.h"

#include <armadillo>

using namespace arma;

Slater::Slater(Config *config, Orbital* orbitals[]) :
    nDimensions(config->nDimensions()),
    nParticles(config->nParticles()),
    orbitals(orbitals)
{
    // the slater determinant can de separated into two parts, one for spin up and the other for spin down
    matrixUpOld = zeros<mat>(nParticles / 2, nParticles / 2);
    matrixDownOld = zeros<mat>(nParticles / 2, nParticles / 2);
}

/*!
  Note: The first half of the particles have spin up, while the others are spin down.
  */
void Slater::constructMatrix(vec2 r[]) {
    for(int i = 0; i < nParticles / 2; i++) {
        for(int j = 0; j < nParticles / 2; j++) {
            matrixUpOld(i,j) = orbitals[j]->evaluate(r[i]);
            matrixDownOld(i,j) = orbitals[j]->evaluate(r[i + nParticles / 2]);
        }
    }
}

double Slater::determinant(vec2 r[]) {
    constructMatrix(r);
    // the slater determinant can be divided into two parts multiplied together
    double detUp = det(matrixUpOld);
    double detDown = det(matrixDownOld);
//    double detUp = 1;
//    double detDown = 2;
    return detUp * detDown;
}

Slater::~Slater() {
}
