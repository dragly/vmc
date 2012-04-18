#ifndef ORBITAL_H
#define ORBITAL_H

#include "../hermite.h"
#include "../config.h"

#include <armadillo>

using namespace arma;

/*!
  A class that contains the single particle wave function.
  */
class Orbital
{
public:
    Orbital(double nx, double ny, Config *config);

    double evaluate(const vec2 r);
    void setParameters(double* parameters);
private:
    // the quantum numbers, one for each dimension in this problem
    double m_nx;
    double m_ny;

    Config *m_config;
    double m_omega;

    double m_alpha;
    double m_beta;
};

#endif // ORBITAL_H
