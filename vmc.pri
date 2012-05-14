myHome = $$system(echo $HOME)
UNAME = $$system(cat /proc/version | grep -o "Red Hat 4.1")
message($$UNAME)
contains(UNAME, Red) {
    message(Old Red Hat Install)
    armadilloDir = $$myHome/apps/armadillo-2.4.4-rhel5
} else {
    message(Newer distribution)
    armadilloDir = $$myHome/apps/armadillo-2.4.4
}
message(Assuming Armadillo is installed under $$armadilloDir/include)
INCLUDEPATH+=$$armadilloDir/include

#LIBS += -llapack
#LIBS += -lblas
#LIBS += -lblas -llapack
LIBS += -llapack -L/usr/lib64/atlas/ -L$$armadilloDir -larmadillo

QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile)
QMAKE_CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
#QMAKE_CXXFLAGS += -std=c++0x

contains(CONFIG,debug) {
    message(Debug mode. Disabling optimization)
    QMAKE_CXXFLAGS -= -O2
    QMAKE_CXXFLAGS += -O0
} else {
    message(Release mode. Optimization at O3.)
    QMAKE_CXXFLAGS -= -O2
    QMAKE_CXXFLAGS += -O3
}
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS

message(compiler: $$QMAKE_CXX)
message(compile flags: $$QMAKE_CXXFLAGS_RELEASE)

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
