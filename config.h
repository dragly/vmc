#ifndef CONFIG_H
#define CONFIG_H

class Config
{
public:
    Config(int rank, int nProcesses, int nDimensions, int nParticles);
    int rank() { return m_rank; }
    int nProcesses() { return m_nProcesses; }
    int nDimensions() { return m_nDimensions; }
    int nParticles() { return m_nParticles; }
private:
    int m_rank;
    int m_nProcesses;
    int m_nDimensions;
    int m_nParticles;
};

#endif // CONFIG_H
