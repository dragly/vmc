#-------------------------------------------------
#
# Project created by QtCreator 2012-02-02T16:18:47
#
#-------------------------------------------------


QT       += core

TARGET = vmc

TEMPLATE = app

include(vmc.pri)

OTHER_FILES += config.ini \
    todo.txt \
    vmc.pri \
    scripts_mytasks.pl

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
} else {
    message(Release mode. Optimization kept at default.)
}
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS

message(compiler: $$QMAKE_CXX)
message(compile flags: $$QMAKE_CXXFLAGS_RELEASE)

# Copy config file to shadow build directory
copyCommand = $(COPY_DIR) $$PWD/config.ini $$OUT_PWD/config.ini
copyFiles.commands = $$copyCommand
copyFiles.target = copyFiles
# update todo list
todoCommand = cd $$PWD; $$PWD/scripts_mytasks.pl > tasks.tasks; cd $$OUT_PWD
todoStuff.commands = $$todoCommand
todoStuff.target = todoStuff
# add dependencies
first.depends = $(first) copyFiles todoStuff
export(first.depends)
export(copyFiles.commands)
export(todoStuff.commands)
QMAKE_EXTRA_TARGETS += first copyFiles todoStuff
