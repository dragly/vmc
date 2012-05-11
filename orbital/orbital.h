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

    double evaluate(const vec2 &r) const;
    void setParameters(double* parameters);
    void gradient(const vec2 &r, vec2 &rGradient) const;
    double laplace(const vec2 &r);
    double nx() { return m_nx; }
    double ny() { return m_ny; }
private:
    // the quantum numbers, one for each dimension in this problem
    double m_nx;
    double m_ny;

    Config *config;
    double omega;

    double alpha;
    double beta;
};

#endif // ORBITAL_H
