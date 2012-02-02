#-------------------------------------------------
#
# Project created by QtCreator 2012-02-02T16:18:47
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = part1
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp

# MPI Settings
QMAKE_CXX = mpic++
QMAKE_CXX_RELEASE = mpic++
QMAKE_CXX_DEBUG = mpic++
QMAKE_LINK = mpic++
QMAKE_CC = mpicc

QMAKE_CFLAGS = $$system(mpicc --showme:compile)
QMAKE_CXXFLAGS = $$system(mpic++ --showme:compile)
QMAKE_LFLAGS = $$system(mpic++ --showme:link)

LIBS += -L/usr/lib/openmpi/
INCLUDEPATH += /usr/include/openmpi/openmpi/ompi/mpi/cxx
