#ifndef ORBITAL_H
#define ORBITAL_H

#include "../hermite.h";
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

    double evaluate(vec2 r);
private:
    // the quantum numbers, one for each dimension in this problem
    double m_nx;
    double m_ny;

    double m_omega;
    Config *m_config;
};

#endif // ORBITAL_H
