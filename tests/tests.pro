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
