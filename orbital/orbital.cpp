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
double Orbital::evaluate(const vec2 &r) const
{
    double sqrtAlphaOmega = sqrt(m_omega * m_alpha);
    double Hx = Hermite::evaluate(m_nx, sqrtAlphaOmega*r[0]);
    double Hy = Hermite::evaluate(m_ny, sqrtAlphaOmega*r[1]);
    return Hx*Hy* exp(-m_alpha*m_omega*(r[0]*r[0] + r[1]*r[1])/2);
}

void Orbital::setParameters(double *parameters)
{
    m_alpha = parameters[0];
    m_beta = parameters[1];
}

void Orbital::gradient(const vec2 &r, vec2 &rGradient) const {
    double sqrtAlphaOmega = sqrt(m_omega * m_alpha);
    double evaluation = exp(-m_alpha*m_omega*(r[0]*r[0] + r[1]*r[1])/2);
    double Hx = Hermite::evaluate(m_nx, sqrtAlphaOmega * r[0]);
    double Hy = Hermite::evaluate(m_ny, sqrtAlphaOmega * r[1]);
    double dHx = Hermite::derivative(m_nx, sqrtAlphaOmega * r[0]);
    double dHy = Hermite::derivative(m_ny, sqrtAlphaOmega * r[1]);
    rGradient[0] = Hy * m_omega * m_alpha * (sqrtAlphaOmega * dHx - Hx * r[0]);
    rGradient[1] = Hx * m_omega * m_alpha * (sqrtAlphaOmega * dHy - Hy * r[1]);
    rGradient[0] *= evaluation;
    rGradient[1] *= evaluation;
}
