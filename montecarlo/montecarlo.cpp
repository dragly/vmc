#include "montecarlo.h"

MonteCarlo::MonteCarlo(WaveFunction *wave, Hamiltonian *hamiltonian) :
    m_wave(wave),
    m_hamiltonian(hamiltonian)
{
}
