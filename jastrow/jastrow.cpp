#include "jastrow.h"

#include "../config.h"

Jastrow::Jastrow(Config *config) :
    m_nParticles(config->nParticles())
{
}

double Jastrow::evaluate(vec2* r) {
    double wf = 1;
    double a = 1; // TODO Should be set to something with regards to spin!
    // TODO Does this need a factor in front of it to account for counting only half of the particles?
    double vec[2];
    for(int i = 0; i < m_nParticles; i++) {
        for(int j = i + 1; j < m_nParticles; j++) {
            vec[0] = r[i][0]-r[j][0];
            vec[1] = r[i][1]-r[j][1];
            double r12 = sqrt(vec[0]*vec[0]+vec[1]*vec[1]);
            double jastrowArgument = (a * r12) / (1 + m_beta * r12);
            wf *= exp(jastrowArgument);
        }
    }
    return wf;
}

void Jastrow::setParameters(double *parameters)
{
    m_beta = parameters[1];
}
