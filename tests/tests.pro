#-------------------------------------------------
#
# Project created by QtCreator 2012-03-08T15:27:59
#
#-------------------------------------------------

QT       += testlib

QT       -= gui

TARGET = tst_vmctests
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

include(../vmc.pri)

MY_SOURCES =

for(source, SOURCES):MY_SOURCES+=../$$source

MY_HEADERS =

for(header, SOURCES):MY_HEADERS+=../$$header

HEADERS = $$MY_HEADERS

SOURCES = $$MY_SOURCES
SOURCES -= ../main.cpp

SOURCES += tst_vmctests.cpp

DEFINES += SRCDIR=\\\"$$PWD/\\\"

LIBS += -larmadillo

QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

#QMAKE_CFLAGS += $$system(mpicc --showme:compile)
#QMAKE_LFLAGS += $$system(mpicxx --showme:link)
#QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile)
#QMAKE_CXXFLAGS += -DMPICH_IGNORE_CXX_SEEK
#QMAKE_LFLAGS += -llapack

#contains(CONFIG,debug) {
#    message(Debug mode. Disabling optimization)
#    QMAKE_CXXFLAGS -= -O2
#    QMAKE_CXXFLAGS += -O0
#} else {
#    message(Release mode. Optimization at O3.)
#    QMAKE_CXXFLAGS -= -O2
#    QMAKE_CXXFLAGS += -O3
#}
#QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
