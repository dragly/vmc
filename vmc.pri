myHome = $$system(echo $HOME)
message(Assuming Armadillo is installed under $$myHome/apps/armadillo-2.4.4/include)
INCLUDEPATH+=$$myHome/apps/armadillo-2.4.4/include

#LIBS += -llapack

HEADERS += \
    matrix.h \
    wavefunction/wavefunction.h \
    montecarlo/montecarlo.h \
    montecarlo/montecarlostandard.h \
    random.h \
    wavefunction/wavesimple.h \
    utils.h \
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
    montecarlo/montecarlometropolishastings.h\
    slater/slater.h \
    minimizer/minimizerevolutionary.h \
    hermite.h \
    orbital/orbital.h \
    minimizer/minimizerstandard.h \
    minimizer/minimizer.h \
    wavefunction/waveslater.h \
    jastrow/jastrow.h

SOURCES += main.cpp \
    matrix.cpp \
    wavefunction/wavefunction.cpp \
    montecarlo/montecarlo.cpp \
    montecarlo/montecarlostandard.cpp \
    random.cpp \
    wavefunction/wavesimple.cpp \
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
    montecarlo/montecarlometropolishastings.cpp\
    slater/slater.cpp \
    minimizer/minimizerevolutionary.cpp \
    hermite.cpp \
    orbital/orbital.cpp \
    minimizer/minimizer.cpp \
    minimizer/minimizerstandard.cpp \
    wavefunction/waveslater.cpp \
    jastrow/jastrow.cpp
