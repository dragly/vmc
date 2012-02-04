#-------------------------------------------------
#
# Project created by QtCreator 2012-02-02T16:18:47
#
#-------------------------------------------------

QT       += core

TARGET = vmc

TEMPLATE = app

SOURCES += main.cpp \
    matrix.cpp \
    wavefunction.cpp \
    wavestandard.cpp \
    montecarlo.cpp \
    montecarlostandard.cpp \
    random.cpp \
    wavesimple.cpp \
    hamiltonian.cpp \
    hamiltonianstandard.cpp \
    hamiltoniansimple.cpp \
    minimizer.cpp \
    minimizerstandard.cpp \
    mainapplication.cpp \
    inih/ini.c \
    inih/cpp/INIReader.cpp

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS = $$system(mpicc --showme:compile)
QMAKE_LFLAGS = $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS = $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
message($$QMAKE_CXXFLAGS_RELEASE)

HEADERS += \
    matrix.h \
    wavefunction.h \
    wavestandard.h \
    montecarlo.h \
    montecarlostandard.h \
    random.h \
    wavesimple.h \
    hamiltonian.h \
    hamiltonianstandard.h \
    utils.h \
    hamiltoniansimple.h \
    minimizer.h \
    minimizerstandard.h \
    mainapplication.h \
    inih/ini.h \
    inih/cpp/INIReader.h

FORMS +=






































