#ifndef CONFIG_H
#define CONFIG_H

#include "wavefunction/wavefunction.h"
#include "hamiltonian/hamiltonian.h"
#include "montecarlo/montecarlo.h"

#include <iostream>

using namespace std;

class Config
{
public:
    Config(int rank, int nProcesses);
    int rank() {
        return m_rank;
    }
    int nProcesses() const  {
        return m_nProcesses;
    }
    int nDimensions() const {
        return m_nDimensions;
    }
    int nParticles() const {
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

    MonteCarlo* monteCarlo() const
    {
        return m_monteCarlo;
    }

    double charge() const
    {
        return m_charge;
    }

    double stepLength() const
    {
        return m_stepLength;
    }

    void loadConfiguration(INIReader *settings);

    void setWave(WaveFunction* arg)
    {
        m_wave = arg;
    }

    void setHamiltonian(Hamiltonian* arg)
    {
        m_hamiltonian = arg;
    }

    void setMonteCarlo(MonteCarlo* arg)
    {
        m_monteCarlo = arg;
    }

    double omega() const
    {
        return m_omega;
    }

    void setOmega(double arg)
    {
        m_omega = arg;
    }

private:
    int m_rank;
    int m_nProcesses;
    int m_nParticles;
    int m_nDimensions;
    double m_charge;
    double m_stepLength;
    WaveFunction* m_wave;
    Hamiltonian* m_hamiltonian;
    MonteCarlo* m_monteCarlo;
    double m_omega;
};

#endif // CONFIG_H
