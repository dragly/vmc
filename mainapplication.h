#ifndef MAINAPPLICATION_H
#define MAINAPPLICATION_H

//#include <QSettings>

class Minimizer;
class INIReader;

class MainApplication
{
public:
    MainApplication(int* argc, char ***argv);

    enum Mode {
        BlockingMode,
        DensityMode,
        MinimizerMode
    };

    void loadConfiguration();
    void runConfiguration();

    void runMinimizer();
    void runBlocking();
    void runDensity();
    void finalize();
private:
    INIReader *settings;



    int* argc;
    char*** argv;
    int rank;
    int nProcesses;

    Mode mode;
};

#endif // MAINAPPLICATION_H
