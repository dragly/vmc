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
    matrix.cpp

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = mpicxx
QMAKE_CXX_DEBUG = mpicxx
QMAKE_LINK = mpicxx
QMAKE_CC = mpicc

QMAKE_CFLAGS = $$system(mpicc --showme:compile)
QMAKE_CXXFLAGS = $$system(mpic++ --showme:compile)
QMAKE_LFLAGS = $$system(mpic++ --showme:link)

HEADERS += \
    matrix.h


