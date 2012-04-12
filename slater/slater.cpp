#include "slater.h"

#include "../config.h"

#include <armadillo>

using namespace arma;

Slater::Slater(Config *config) :
    m_nDimensions(config->nDimensions()),
    m_nParticles(config->nParticles())
{
    matrix = new mat(m_nParticles / 2, m_nParticles / 2);
}

Slater::~Slater() {
    delete matrix;
}
