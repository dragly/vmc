#include "slater.h"

#include "../config.h"
#include "../orbital/orbital.h"

#include <armadillo>

using namespace arma;

Slater::Slater(Config *config) :
    m_nDimensions(config->nDimensions()),
    m_nParticles(config->nParticles())
{
    // the slater determinant can de separated into two parts, one for spin up and the other for spin down
    matrixUp = zeros<mat>(m_nParticles / 2, m_nParticles / 2);
    matrixDown = zeros<mat>(m_nParticles / 2, m_nParticles / 2);
}

void Slater::constructMatrix(vec2 *&r, Orbital* orbitals[]) {
    for(int i = 0; i < m_nParticles / 2; i++) {
        for(int j = 0; j < m_nParticles / 2; j++) {
            matrixUp(i,j) = orbitals[j]->evaluate(r[i]);
            matrixDown(i,j) = orbitals[j]->evaluate(r[i + m_nParticles / 2]);
        }
    }
}

double Slater::determinant(vec2 *&r, Orbital* orbitals[]) {
    constructMatrix(r, orbitals);
    // the slater determinant can be divided into two parts multiplied together
    double detUp = det(matrixUp);
    double detDown = det(matrixDown);
    double normFactorial = 1;
    for(int i = 2; i < m_nParticles / 2; i++) {
        normFactorial *= i;
    }
    return 1 / sqrt(normFactorial) * detUp * detDown;
}

Slater::~Slater() {
}
