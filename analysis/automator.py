#!/usr/bin/python
# -*- coding: utf-8 -*-
import subprocess
import dircache
import inspect, os
import ConfigParser
from sys import argv

from automatorsettings import *

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

autoMode = argv[1]
runLocal = False

currentDirectory = os.path.split(os.path.abspath(__file__))[0]
if len(argv) > 2:
    runLocal = argv[2]

servers = servers.split(" ")

files = dircache.listdir(runFolder)
i = 0

if autoMode == "kill":
    for server in servers:
	print "Calling: ssh " + server +  " 'killall vmc'"
	subprocess.call("ssh " + server +  " 'killall vmc'", shell=True)
elif autoMode == "status":
    for server in servers:
	print server + ":"
	nProcs = subprocess.check_output("ssh " + server +  " 'ps aux | grep generator.py | wc -l'", shell=True)
	if(int(nProcs) > 2):
	    print bcolors.OKGREEN + "Running " + nProcs.replace("\n","") + bcolors.ENDC
	users = subprocess.check_output("ssh " + server +  " 'users'", shell=True)
	if users != "":
	    print bcolors.WARNING + users.replace("\n","") + bcolors.ENDC
elif autoMode == "run":
    for afile in files:
	print "Found job in " + afile
	configName = afile
	
	myDir = runFolder + "/" + configName + "/"
	
	config = ConfigParser.ConfigParser()
	configFileName = myDir + "config.ini"
	config.read(configFileName)
	mode = config.get("General", "mode")
	fileop = open(configFileName)
	configData = fileop.read()
	
	prevConfigFile = myDir + prevrunConfigFileName
	prevConfigData = ""
	if os.path.exists(prevConfigFile):
	    prevConfigFileOp = open(prevConfigFile)
	    prevConfigData = prevConfigFileOp.read()
	    
	if prevConfigData == configData:
	    print bcolors.WARNING + "Config is the same. No need to run again." + bcolors.ENDC
	elif os.path.exists(myDir + lockFileName):
	    lockFileOp = open(myDir + lockFileName)
	    lockHost = lockFileOp.read()
	    lockFileOp.close()
	    print bcolors.WARNING + "Run is already locked by process on host " + lockHost + bcolors.ENDC
	else:
	    runCommand = "python generator.py " + configName
	    serverOccupied = True
	    while(serverOccupied):
		nProcs = subprocess.check_output("ssh " + servers[i] +  " 'ps aux | grep generator.py | wc -l'", shell=True)
		if(int(nProcs) < 3):
		    serverOccupied = False
		else:
		    print bcolors.WARNING + "Server " + servers[i] + " already in use. Choosing next." + bcolors.ENDC
		    i += 1
		    if i >= len(servers):
			print "All servers occupied. Exiting!"
			exit(0)
			break;
	    if runLocal:
		print "Running " + configName + " locally"
		subprocess.call(runCommand, shell=True)
	    else:
		print bcolors.OKGREEN + "Launching job for config " + configName + " on " + servers[i] + bcolors.ENDC
		subprocess.Popen("ssh " + servers[i] +  " 'cd "+ currentDirectory + " && " + runCommand + " &> output-" + configName + " && echo " + servers[i] + " done'", shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
		i += 1
            print "You may watch the output from the servers in the output-* files."
elif autoMode == "plot":
    minimizerList = ""
    densityList = ""
    for afile in files:
	configName = afile
	myDir = runFolder + "/" + configName + "/"
	config = ConfigParser.ConfigParser()
	configFileName = myDir + configFileName
	config.read(configFileName)
	mode = config.get("General", "mode")
	if mode == "minimizer":
	    minimizerList += myDir + " "
	elif mode == "density":
	    densityList += myDir + " "
	else:
	    print "Unknown mode"
	    
    subprocess.call("python minimizerplot.py " + minimizerList, shell=True)
    subprocess.call("python densityplot.py " + densityList, shell=True)
print "Automator finished."
