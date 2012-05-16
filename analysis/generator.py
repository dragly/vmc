#!/usr/bin/python
# -*- coding: utf-8 -*-
from sys import argv
import inspect, os
import time 
import ConfigParser
import subprocess

from automatorsettings import *

configName = argv[1]

myDir = runFolder + "/" + configName + "/"

lockFile = open(myDir + lockFileName,"w+")
lockFile.write(os.uname()[1])
lockFile.close()

config = ConfigParser.ConfigParser()
config.readfp(open(runFolder + "/" + configName + "/" + configFileName))
mode = config.get("General", "mode")
if mode == "minimizer" or mode == "density":
    runCommand = runCommandMPI   
nCycles = config.getint("MinimizerStandard", "nCycles")

if not os.path.exists(myDir):
    os.makedirs(myDir)
    
print "Run started by " + os.uname()[1] + " at " + str(time.time()) + " with " + str(nCycles) + " cycles\n"

subprocess.call("./buildandrun.sh " + myDir + " " + runCommand + " " + configFileName + " " + prevrunConfigFileName, shell=True)

# this script is not included, but basically just calls a command that sends me an e-mail when the job is done. May be commented out
#subprocess.call("./sendmail.script", shell=True)

os.remove(myDir + "lock")
print "Run done at " + time.strftime("%Y-%m-%d %H:%M:%S")