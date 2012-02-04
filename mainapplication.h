#ifndef MAINAPPLICATION_H
#define MAINAPPLICATION_H

#include <QSettings>

class Minimizer;

class MainApplication
{
public:
    MainApplication(int* argc, char ***argv);

    void loadConfiguration();
    int runApplication();

private:
    QSettings *settings;
    Minimizer *minimizer;

    int* argc;
    char*** argv;
    int rank;
    int nProcesses;
};

#endif // MAINAPPLICATION_H
