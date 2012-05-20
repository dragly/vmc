#ifndef CONFIG_H
#define CONFIG_H

#include "wavefunction/wavefunction.h"
#include "hamiltonian/hamiltonian.h"
#include "montecarlo/montecarlo.h"

#include <iostream>


class INIParser;
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
    void setNDimensions(int n) {
        m_nDimensions = n;
    }

    int nParticles() const {
        return m_nParticles;
    }

    void setNParticles(int n) {
        m_nParticles = n;
    }

    WaveFunction* wave() const
    {
        return m_wave;
    }

    Hamiltonian* hamiltonian() const
    {
        return m_hamiltonian;
    }

    string monteCarloClass() const
    {
        return m_monteCarloClass;
    }

    double stepLength() const
    {
        return m_stepLength;
    }

    void loadConfiguration(INIParser *settings);

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

    void setTau(double arg) {
        m_tau = arg;
    }

    double tau() {
        return m_tau;
    }

    void setDiffusionConstant(double arg) {
        m_diffusionConstant = arg;
    }

    double diffusionConstant() {
        return m_diffusionConstant;
    }

    bool interactionEnabled() const
    {
        return m_interactionEnabled;
    }

    void setInteractionEnabled(bool arg)
    {
        m_interactionEnabled = arg;
    }

    long *idum() {
        return &m_idum;
    }

private:
    int m_rank;
    int m_nProcesses;
    int m_nParticles;
    int m_nDimensions;
    double m_stepLength;
    MonteCarlo* m_monteCarlo;
    WaveFunction* m_wave;
    Hamiltonian* m_hamiltonian;
    string m_monteCarloClass;
    double m_omega;
    bool m_interactionEnabled;
    long m_idum;
    double m_tau;
    double m_diffusionConstant;
};

#endif // CONFIG_H
