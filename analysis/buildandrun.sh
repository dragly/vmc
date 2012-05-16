#!/bin/bash
# The following input parameters must be set as arguments, i.e. ./buildandrun.sh $1 $2 $3 $4
# $1 = folder to run in
# $2 = run command
# $3 = configFileName
# $4 = prevrunConfigFileName
echo Building on $HOSTNAME;
cd $1;
echo Building in $(pwd)
mpd &
PATH=$PATH:/usr/lib64/qt-3.3/bin/ qmake CONFIG-=debug CONFIG+=nocopy ../../../vmc.pro && make clean && make && LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/apps/armadillo-2.4.4-rhel5/ $2 && cp $3 $4
