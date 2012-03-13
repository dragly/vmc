#include "config.h"

Config::Config(int rank, int nProcesses, int nDimensions, int nParticles, double charge, double stepLength, WaveFunction *wave, Hamiltonian *hamiltonian) :
    m_rank(rank),
    m_nProcesses(nProcesses),
    m_nDimensions(nDimensions),
    m_nParticles(nParticles)

{
}
