#ifndef ONERUN_H
#define ONERUN_H

class INIParser;
class WaveFunction;
class MonteCarlo;
class Hamiltonian;
class Config;
class Blocker;

#include <string>

/*!
  * \brief Runs a variational Monte Carlo run with a given MonteCarlo class for a set of parameters
  */
class OneRun
{
public:
    OneRun(Config *config);

    void loadConfiguration(INIParser *settings);
    void run();
    void writeBlockData();
    ~OneRun();
private:
    double nSamples;
    double alpha;
    double beta;
    int myRank;
    int m_nProcesses;
    double *allEnergies;
    WaveFunction *wave;
    MonteCarlo *monteCarlo;
    Hamiltonian *hamiltonian;
    Config *config;
    Blocker *blocker;
    bool onlyBlocking;
    std::string scratchDir;
};

#endif // ONERUN_H
