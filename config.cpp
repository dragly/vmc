#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

//#include "inih/cpp/INIReader.h"
#include "inih/ini.h"
#include "wavefunction/waveslater.h"
#include "hamiltonian/hamiltonianideal.h"
#include "montecarlo/montecarlostandard.h"

using namespace std;

Config::Config(int rank, int nProcesses) :
    m_rank(rank),
    m_nProcesses(nProcesses),
    m_nParticles(2),
    m_nDimensions(2),
    m_stepLength(1.0),
    m_wave(0),
    m_hamiltonian(0),
    m_monteCarloClass("MonteCarloStandard"),
    m_omega(1),
    m_interactionEnabled(true),
    m_idum(-rank*time(NULL)),
    m_diffusionConstant(0.5),
    m_tau(0.01)
{
}

void Config::loadConfiguration(INIParser* settings) {
    m_nParticles = settings->GetInteger("General", "nParticles", m_nParticles);
    m_nDimensions = settings->GetInteger("General", "nDimensions", m_nDimensions);
    m_interactionEnabled = settings->GetBoolean("General", "interactionEnabled", m_interactionEnabled);
    m_stepLength = atof(settings->Get("General", "stepLength", "1.0").c_str());
    m_omega = atof(settings->Get("General", "omega", "1.0").c_str());

    // Wave properties
    string waveClass = settings->Get("Wave","class", "WaveSimple");
    std::cout << waveClass << std::endl;
    m_wave = WaveFunction::fromName(waveClass, this);
//    m_wave = new WaveSlater(this);
    if(m_wave == 0) {
        cerr << "Unknown wave class '" << waveClass << "'" << endl;
        exit(99);
    }
    m_wave->loadConfiguration(settings);

    // Hamiltonian
    string hamiltonianClass = settings->Get("Hamiltonian","class", "HamiltonianSimple");
    std::cout << hamiltonianClass << std::endl;
    m_hamiltonian = Hamiltonian::fromName(hamiltonianClass, this);
//    m_hamiltonian = new HamiltonianIdeal(this);
    if(m_hamiltonian == 0) {
        cerr << "Unknown hamiltonian class '" << hamiltonianClass << "'" << endl;
        exit(98);
    }

    // Monte Carlo sampler
    m_monteCarloClass = settings->Get("MonteCarlo","class", "MonteCarloStandard");
    std::cout << m_monteCarloClass << std::endl;
    m_monteCarlo = MonteCarlo::fromName(m_monteCarloClass, this);
//    m_monteCarlo = new MonteCarloStandard(this);
    if(m_monteCarlo == 0) {
        cerr << "Unknown Monte Carlo class '" << m_monteCarloClass << "'" << endl;
        exit(97);
    }
}
