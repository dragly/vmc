#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include "inih/cpp/INIReader.h"

using namespace std;

Config::Config(int rank, int nProcesses) :
    m_rank(rank),
    m_nProcesses(nProcesses),
    m_nParticles(2),
    m_nDimensions(2),
    m_charge(1.0),
    m_stepLength(1.0),
    m_wave(0),
    m_hamiltonian(0),
    m_monteCarlo(0)
{
}

void Config::loadConfiguration(INIReader* settings) {
    m_nParticles = settings->GetInteger("General", "nParticles", 2);
    m_nDimensions = settings->GetInteger("General", "nDimensions", 2);
    m_charge = atof(settings->Get("General", "charge", "1.0").c_str());
    m_stepLength = atof(settings->Get("General", "stepLength", "1.0").c_str());

    // Wave properties
    string waveClass = settings->Get("Wave","class", "WaveSimple");
    m_wave = WaveFunction::fromName(waveClass, this);
    if(m_wave == 0) {
        cerr << "Unknown wave class '" << waveClass << "'" << endl;
        exit(99);
    }
    m_wave->loadConfiguration(settings);

    // Hamiltonian
    string hamiltonianClass = settings->Get("Hamiltonian","class", "HamiltonianSimple");
    m_hamiltonian = Hamiltonian::fromName(hamiltonianClass, this, m_charge);
    if(m_hamiltonian == 0) {
        cerr << "Unknown hamiltonian class '" << hamiltonianClass << "'" << endl;
        exit(98);
    }

    // Monte Carlo sampler
    string monteCarloClass = settings->Get("MonteCarlo","class", "MonteCarloSimple");
    m_monteCarlo = MonteCarlo::fromName(monteCarloClass, this);
    if(m_monteCarlo == 0) {
        cerr << "Unknown Monte Carlo class '" << monteCarloClass << "'" << endl;
        exit(98);
    }
}
