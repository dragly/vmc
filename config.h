#ifndef CONFIG_H
#define CONFIG_H

#include "wavefunction.h"
#include "hamiltonian/hamiltonian.h"

class Config
{
    Q_PROPERTY(double stepLength READ stepLength)
public:
    Config(int rank, int nProcesses, int nDimensions, int nParticles, double charge, double stepLength, WaveFunction* wave, Hamiltonian* hamiltonian);
    int rank() {
        return m_rank;
    }
    int nProcesses() {
        return m_nProcesses;
    }
    int nDimensions() {
        return m_nDimensions;
    }
    int nParticles() {
        return m_nParticles;
    }

    WaveFunction* wave() const
    {
        return m_wave;
    }

    Hamiltonian* hamiltonian() const
    {
        return m_hamiltonian;
    }

    double charge() const
    {
        return m_charge;
    }

    double stepLength() const
    {
        return m_stepLength;
    }

private:
    int m_rank;
    int m_nProcesses;
    int m_nDimensions;
    int m_nParticles;
    WaveFunction* m_wave;
    Hamiltonian* m_hamiltonian;
    double m_charge;
    double m_stepLength;
};

#endif // CONFIG_H
