myHome = $$system(echo $HOME)
INCLUDEPATH+=$$myHome/apps/armadillo-2.4.4/include

SOURCES += main.cpp \
    matrix.cpp \
    wavefunction/wavefunction.cpp \
    montecarlo/montecarlo.cpp \
    montecarlo/montecarlostandard.cpp \
    random.cpp \
    wavefunction/wavesimple.cpp \
    minimizer.cpp \
    minimizerstandard.cpp \
    mainapplication.cpp \
    inih/ini.c \
    inih/cpp/INIReader.cpp \
    wavefunction/waveideal.cpp \
    hamiltonian/hamiltonianideal.cpp \
    hamiltonian/hamiltoniansimple.cpp \
    hamiltonian/hamiltonian.cpp \
    blocker.cpp \
    densityplotter.cpp \
    config.cpp \
    montecarlo/montecarlometropolishastings.cpp

HEADERS += \
    matrix.h \
    wavefunction/wavefunction.h \
    wavefunction/wavestandard.h \
    montecarlo/montecarlo.h \
    montecarlo/montecarlostandard.h \
    random.h \
    wavefunction/wavesimple.h \
    utils.h \
    minimizer.h \
    minimizerstandard.h \
    mainapplication.h \
    inih/ini.h \
    inih/cpp/INIReader.h \
    wavefunction/waveideal.h \
    hamiltonian/hamiltonianideal.h \
    hamiltonian/hamiltonian.h \
    hamiltonian/hamiltoniansimple.h \
    hamiltonian/hamiltonianstandard.h \
    blocker.h \
    densityplotter.h \
    config.h \
    montecarlo/montecarlometropolishastings.h

