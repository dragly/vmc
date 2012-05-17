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
    // Optimizing calculations that are much reused (23 % self cost decrease)
    double alphaOmega = omega * alpha;
    double sqrtAlphaOmega = sqrt(alphaOmega);
    double sqrtAlphaOmegax = sqrtAlphaOmega * x;
    double sqrtAlphaOmegay = sqrtAlphaOmega * y;
    double alphaOmegaalphaOmega = alphaOmega * alphaOmega;
    // All calculations
    double evaluation = exp(-alphaOmega*(x*x + y*y)/2);
    double Hx = Hermite::evaluate(m_nx, sqrtAlphaOmegax);
    double Hy = Hermite::evaluate(m_ny, sqrtAlphaOmegay);
    double dHx = Hermite::derivative(m_nx, sqrtAlphaOmegax);
    double dHy = Hermite::derivative(m_ny, sqrtAlphaOmegay);
    double ddHx = Hermite::doubleDerivative(m_nx, sqrtAlphaOmegax);
    double ddHy = Hermite::doubleDerivative(m_ny, sqrtAlphaOmegay);
    double doubleDerivativeX = Hy * evaluation * ( alphaOmega * ddHx - sqrtAlphaOmegax * dHx * alphaOmega
                                                   - sqrtAlphaOmegax * dHx * alphaOmega
                                                   - Hx * alphaOmega + Hx * alphaOmegaalphaOmega * x * x);
//    std::cout << doubleDerivativeX << std::endl;
    double doubleDerivativeY = Hx * evaluation * ( alphaOmega * ddHy - sqrtAlphaOmegay * dHy * alphaOmega
                                                   - sqrtAlphaOmegay * dHy * alphaOmega
                                                   - Hy * alphaOmega + Hy * alphaOmegaalphaOmega * y * y);
//    std::cout << doubleDerivativeY << std::endl;
    return doubleDerivativeX + doubleDerivativeY;
}
