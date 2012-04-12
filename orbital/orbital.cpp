#include "orbital.h"

Orbital::Orbital(double nx, double ny, Config *config) :
    m_nx(nx),
    m_ny(ny),
    m_config(config),
    m_omega(config->omega()),
    m_alpha(1),
    m_beta(1)
{
}

/*!
  Evaluates the orbital at the given position.

  @param r Position to evaluate
  @returns Value of the orbital at the given position
  */
double Orbital::evaluate(vec2 r)
{
    double sqrtomega = sqrt(m_omega);
    double Hx = Hermite::evaluate(r[0], sqrtomega*m_nx);
    double Hy = Hermite::evaluate(r[0], sqrtomega*m_ny);
    return Hx*Hy* exp(-m_alpha*m_omega*(r[0]*r[0] + r[1]*r[1])/2);
}

void Orbital::setParameters(double alpha, double beta)
{
    m_alpha = alpha;
    m_beta = beta;
}

