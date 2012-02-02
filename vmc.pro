#-------------------------------------------------
#
# Project created by QtCreator 2012-02-02T16:18:47
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = vmc
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp \
    matrix.cpp \
    wavefunction.cpp \
    wavestandard.cpp \
    montecarlo.cpp \
    montecarlostandard.cpp \
    random.cpp

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = mpicxx
QMAKE_CXX_DEBUG = mpicxx
QMAKE_LINK = mpicxx
QMAKE_CC = mpicc

QMAKE_CFLAGS = $$system(mpicc --showme:compile)
QMAKE_LFLAGS = $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS = $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE = $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
message($$QMAKE_CXXFLAGS_RELEASE)

HEADERS += \
    matrix.h \
    wavefunction.h \
    wavestandard.h \
    montecarlo.h \
    montecarlostandard.h \
    random.h












