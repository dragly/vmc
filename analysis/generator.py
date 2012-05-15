#!/usr/bin/python
# -*- coding: utf-8 -*-
from sys import argv
import inspect, os
import time 
import ConfigParser
import subprocess

configName = argv[1]

myDir = "runs/" + configName + "/"

lockFile = open(myDir + "lock","w+")
lockFile.write(os.uname()[1])
lockFile.close()

config = ConfigParser.ConfigParser()
config.readfp(open('runs/' + configName + "/config.ini"))
mode = config.get("General", "mode")
nCycles = config.getint("MinimizerStandard", "nCycles")
if not os.path.exists(myDir):
    os.makedirs(myDir)
print "Run started by " + os.uname()[1] + " at " + str(time.time()) + " with " + str(nCycles) + " cycles\n"
subprocess.call("echo I am $HOSTNAME", shell=True)
commandPrefix = ""
if mode == "minimizer" or mode == "density":
    commandPrefix = "'mpirun -n 4'"
subprocess.call("./buildandrun.sh " + configName + " " + commandPrefix, shell=True)

# this script is not included, but basically just calls a command that sends me an e-mail when the job is done. May be commented out
#subprocess.call("./sendmail.script", shell=True)
os.remove(myDir + "lock")
print "Done " + time.strftime("%Y-%m-%d %H:%M:%S")