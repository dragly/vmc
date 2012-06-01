#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

//#include "inih/cpp/INIReader.h"
#include "inih/ini.h"
#include "wavefunction/waveslater.h"
#include "hamiltonian/hamiltonianideal.h"
#include "montecarlo/standardmontecarlo.h"

using namespace std;

Config::Config(int myRank, int m_nProcesses) :
    m_rank(myRank),
    m_m_nProcesses(m_nProcesses),
    m_nParticles(2),
    m_nDimensions(2),
    m_stepLength(1.0),
    m_wave(0),
    m_hamiltonian(0),
    m_monteCarloClass("MonteCarloStandard"),
    m_omega(1.0),
    m_interactionEnabled(true),
    m_idum(-(1 + myRank)*time(NULL) / 100000), // TODO figure out why this must be divided by 1000000
    m_diffusionConstant(0.5)
{
}

void Config::loadConfiguration(INIParser* settings) {
    m_nParticles = settings->GetDouble("General", "nParticles", m_nParticles);
    m_nDimensions = settings->GetDouble("General", "nDimensions", m_nDimensions);
    m_interactionEnabled = settings->GetBoolean("General", "interactionEnabled", m_interactionEnabled);
    m_stepLength = settings->GetDouble("General", "stepLength", m_stepLength);
    m_omega = settings->GetDouble("General", "omega", m_omega);

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
    if(m_monteCarloClass == "MonteCarloMetropolisHastings" && m_stepLength > 0.05) {
        std::cerr << "WARNING! You are running importance sampling with a big step length. This will likely cause errors in your results." << std::endl;
    }
//    m_monteCarlo = new MonteCarloStandard(this);
    if(m_monteCarlo == 0) {
        cerr << "Unknown Monte Carlo class '" << m_monteCarloClass << "'" << endl;
        exit(97);
    }
    m_monteCarlo->loadConfiguration(settings);
}
