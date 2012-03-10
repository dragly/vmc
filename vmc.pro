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
    montecarlo.cpp \
    montecarlostandard.cpp \
    random.cpp \
    wavesimple.cpp \
    minimizer.cpp \
    minimizerstandard.cpp \
    mainapplication.cpp \
    inih/ini.c \
    inih/cpp/INIReader.cpp \
    waveideal.cpp \
    hamiltonian/hamiltonianideal.cpp \
    hamiltonian/hamiltoniansimple.cpp \
    hamiltonian/hamiltonian.cpp \
    blocker.cpp \
    densityplotter.cpp \
    config.cpp

OTHER_FILES += config.ini \
    todo.txt \
    vmc.pri

HEADERS += \
    matrix.h \
    wavefunction.h \
    wavestandard.h \
    montecarlo.h \
    montecarlostandard.h \
    random.h \
    wavesimple.h \
    utils.h \
    minimizer.h \
    minimizerstandard.h \
    mainapplication.h \
    inih/ini.h \
    inih/cpp/INIReader.h \
    waveideal.h \
    hamiltonian/hamiltonianideal.h \
    hamiltonian/hamiltonian.h \
    hamiltonian/hamiltoniansimple.h \
    hamiltonian/hamiltonianstandard.h \
    blocker.h \
    densityplotter.h \
    config.h

FORMS +=

# MPI Settings

contains(CONFIG,mpi) {
    message(Using MPI)
    QMAKE_CXX = mpicxx
    QMAKE_CXX_RELEASE = $$QMAKE_CXX
    QMAKE_CXX_DEBUG = $$QMAKE_CXX
    QMAKE_LINK = $$QMAKE_CXX
    QMAKE_CC = mpicc

#    QMAKE_CFLAGS = $$system(mpicc --showme:compile)
#    QMAKE_LFLAGS = $$system(mpicxx --showme:link)
#    QMAKE_CXXFLAGS = $$system(mpicxx --showme:compile)
    QMAKE_CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK

    DEFINES += USE_MPI
} else {
    message(Not using MPI)
    # slow debug mode
}

contains(CONFIG,debug) {
    message(Debug mode. Disabling optimization)
    QMAKE_CXXFLAGS = -O0
}
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS

message(compiler: $$QMAKE_CXX)
message(compile flags: $$QMAKE_CXXFLAGS_RELEASE)

# Copy config file to shadow build directory
copyCommand = $(COPY_DIR) $$PWD/config.ini $$OUT_PWD/config.ini
copyFiles.commands = $$copyCommand
copyFiles.target = copyFiles
first.depends = $(first) copyFiles
export(first.depends)
export(copyFiles.commands)
QMAKE_EXTRA_TARGETS += first copyFiles






































