#include "config.h"

Config::Config(int rank, int nProcesses, int nDimensions, int nParticles) :
    m_rank(rank),
    m_nProcesses(nProcesses),
    m_nDimensions(nDimensions),
    m_nParticles(nParticles)
{
}
