#ifndef MAINAPPLICATION_H
#define MAINAPPLICATION_H

//#include <QSettings>

class Minimizer;
class INIParser;
class Config;

class MainApplication
{
public:
    MainApplication(int* argc, char ***argv);

    enum Mode {
        BlockingMode,
        DensityMode,
        MinimizerMode,
        GeneticMode,
        DiffusionMode,
        OneRunMode
    };

    void loadConfiguration();
    void runConfiguration();

    void finalize();
private:
    INIParser *m_settings;



    int* argc;
    char*** argv;
    int m_rank;
    int m_nProcesses;
    int m_nParticles;
    int m_nDimensions;
    Config *m_config;

    Mode m_mode;
};

#endif // MAINAPPLICATION_H
