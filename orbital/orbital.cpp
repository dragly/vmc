#include "orbital.h"

Orbital::Orbital(double nx, double ny, Config *config) :
    m_nx(nx),
    m_ny(ny),
    config(config),
    omega(config->omega()),
    alpha(1),
    beta(1)
{
}

/*!
  Evaluates the orbital at the given position.

  @param r Position to evaluate
  @returns Value of the orbital at the given position
  */
double Orbital::evaluate(const vec2 &r) const
{
    double sqrtAlphaOmega = sqrt(omega * alpha);
    double Hx = Hermite::evaluate(m_nx, sqrtAlphaOmega*r[0]);
    double Hy = Hermite::evaluate(m_ny, sqrtAlphaOmega*r[1]);
    return Hx*Hy* exp(-alpha*omega*(r[0]*r[0] + r[1]*r[1])/2);
}

void Orbital::setParameters(double *parameters)
{
    alpha = parameters[0];
    beta = parameters[1];
}

void Orbital::gradient(const vec2 &r, vec2 &rGradient) const {
    double x = r[0];
    double y = r[1];
    double sqrtAlphaOmega = sqrt(omega * alpha);
    double evaluation = exp(-alpha*omega*(x*x + y*y)/2);
    double Hx = Hermite::evaluate(m_nx, sqrtAlphaOmega * x);
    double Hy = Hermite::evaluate(m_ny, sqrtAlphaOmega * y);
    double dHx = Hermite::derivative(m_nx, sqrtAlphaOmega * x);
    double dHy = Hermite::derivative(m_ny, sqrtAlphaOmega * y);
    rGradient[0] = Hy * (sqrtAlphaOmega * dHx - Hx * omega * alpha * x);
    rGradient[1] = Hx * (sqrtAlphaOmega * dHy - Hy * omega * alpha * y);
    rGradient[0] *= evaluation;
    rGradient[1] *= evaluation;
}

double Orbital::laplace(const vec2 &r) {
    double x = r[0];
    double y = r[1];
    double sqrtAlphaOmega = sqrt(omega * alpha);
    double evaluation = exp(-alpha*omega*(x*x + y*y)/2);
    double Hx = Hermite::evaluate(m_nx, sqrtAlphaOmega * x);
    double Hy = Hermite::evaluate(m_ny, sqrtAlphaOmega * y);
    double dHx = Hermite::derivative(m_nx, sqrtAlphaOmega * x);
    double dHy = Hermite::derivative(m_ny, sqrtAlphaOmega * y);
    double ddHx = Hermite::doubleDerivative(m_nx, sqrtAlphaOmega * x);
    double ddHy = Hermite::doubleDerivative(m_ny, sqrtAlphaOmega * y);
    double doubleDerivativeX = Hy * evaluation * ( alpha * omega * ddHx - sqrtAlphaOmega * dHx * alpha * omega * x
                                                   - sqrtAlphaOmega * dHx * alpha * omega * x
                                                   - Hx * alpha*omega + Hx * alpha * alpha * omega * omega * x * x);
//    std::cout << doubleDerivativeX << std::endl;
    double doubleDerivativeY = Hx * evaluation * ( alpha * omega * ddHy - sqrtAlphaOmega * dHy * alpha * omega * y
                                                   - sqrtAlphaOmega * dHy * alpha * omega * y
                                                   - Hy * alpha*omega + Hy * alpha * alpha * omega * omega * y * y);
//    std::cout << doubleDerivativeY << std::endl;
    return doubleDerivativeX + doubleDerivativeY;
}
