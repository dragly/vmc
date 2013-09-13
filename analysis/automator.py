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

#    runLocal = argv[2]

servers = servers.split(" ")

files = dircache.listdir(runFolder)
if len(argv) > 2:
    files = argv[2:]
i = 0

if autoMode == "kill":
    for server in servers:
	print "Calling: ssh " + server +  " 'killall vmc'"
	subprocess.call("ssh " + server +  " 'killall vmc'", shell=True)
elif autoMode == "status":
    for afile in files:
        configName = afile
        
        myDir = runFolder + "/" + configName + "/"
        try:
            config = ConfigParser.ConfigParser()
            configFileName = myDir + "config.ini"
            config.read(configFileName)
        except IOError as e:
            print bcolors.FAIL + "Could not find configfile for " + myDir + ". Skipping..."
            continue
        mode = config.get("General", "mode")
        fileop = open(configFileName)
        configData = fileop.read()
        
        prevConfigFile = myDir + prevrunConfigFileName
        prevConfigData = ""
        if os.path.exists(prevConfigFile):
            prevConfigFileOp = open(prevConfigFile)
            prevConfigData = prevConfigFileOp.read()
            
        if prevConfigData == configData:
            #print bcolors.WARNING + "Config is the same. No need to run again." + bcolors.ENDC
            aasfasf = 0
        elif os.path.exists(myDir + lockFileName):
            print "Found job in " + afile
            lockFileOp = open(myDir + lockFileName)
            lockHost = lockFileOp.read()
            lockFileOp.close()
            print bcolors.WARNING + "Run is locked by process on host " + lockHost + bcolors.ENDC
        else:
            print "Found job in " + afile
            print bcolors.OKGREEN + "Updated config file. It needs a rerun." + bcolors.ENDC
            
    for server in servers:
        print server + ":"
        nProcs = subprocess.check_output("ssh " + server +  " 'ps aux | grep vmc | grep svennard | wc -l'", shell=True)
        if(int(nProcs) > 2):
            print bcolors.OKGREEN + "Running " + nProcs.replace("\n","") + bcolors.ENDC
        runstring = "ssh " + server +  " \"ps aux | grep -v 'svennard\|root\|dbus\|rpc\|dbus\|hald\|xfs\|avahi\|gdm\|USER\|ntp' |  awk '{print \$1,\$3}'\""
        users = subprocess.check_output(runstring, shell=True)
        userlist = users.split("\n")
        useroutput = ""
        for userdata in userlist:
            usersplit = userdata.split(" ")
            if len(usersplit) > 1:
                username = usersplit[0]
                cpuusage = float(usersplit[1])
                if cpuusage > 10:
                    useroutput += username + " " + str(cpuusage) + " "
        if useroutput != "":
            print bcolors.FAIL + useroutput + bcolors.ENDC
        #if users != "":
        #    print bcolors.WARNING + users.replace("\n"," ") + bcolors.ENDC
        #users2 = subprocess.check_output("ssh " + server + " \"users | grep -v svennard\"", shell=True)
        #if users2 != "":
        #    print bcolors.FAIL + users2.replace("\n","") + bcolors.ENDC
elif autoMode == "run":
    for afile in files:
	configName = afile
	if len(argv) > 2:
	    configName = afile.replace("runs/", "")
	    configName = configName.replace("/", "")
	else:
	    configName = afile
	    
	myDir = runFolder + "/" + configName + "/"
	
	print "Found job in " + myDir
	    
	try:
            config = ConfigParser.ConfigParser()
            configFileName = myDir + "config.ini"
            config.read(configFileName)
        except IOError as e:
            print bcolors.FAIL + "Could not find configfile for " + myDir + ". Skipping..."
            continue
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
                if i >= len(servers):
                    print "All servers occupied. Exiting!"
                    exit(0)
                    break;
		nProcs = subprocess.check_output("ssh " + servers[i] +  " 'ps aux | grep generator.py | wc -l'", shell=True)
		if(int(nProcs) < 3):
		    serverOccupied = False
		else:
		    print bcolors.WARNING + "Server " + servers[i] + " already in use. Choosing next." + bcolors.ENDC
		    i += 1
	    if runLocal:
		print "Running " + configName + " locally"
		subprocess.call(runCommand, shell=True)
	    else:
		print bcolors.OKGREEN + "Launching job for config " + configName + " on " + servers[i] + bcolors.ENDC
		subprocess.Popen("ssh " + servers[i] +  " 'cd "+ currentDirectory + " && " + runCommand + " &> output/output-" + configName + " && echo " + servers[i] + " done with " + configName + "'", shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
		i += 1
            print "You may watch the output from this job in the output-" + configName + " file."
elif autoMode == "plot":
    plotAll = False
    if len(argv) < 3 or argv[2] == "all":
        plotAll = True
    minimizerList = ""
    densityList = ""
    geneticList = ""
    onerunList = ""
    diffusionList = ""
    distributionList = ""
    for afile in files:
        configName = afile
        myDir = runFolder + "/" + configName + "/"
        if os.path.exists(myDir + lockFileName):
            print "Found job in " + afile
            lockFileOp = open(myDir + lockFileName)
            lockHost = lockFileOp.read()
            lockFileOp.close()
            print bcolors.WARNING + "Run is locked by process on host " + lockHost + ". Skipping locked plots." + bcolors.ENDC
        else:
            config = ConfigParser.ConfigParser()
            configFilePath = myDir + configFileName
            config.read(configFilePath)
            mode = config.get("General", "mode")
            if mode == "minimizer":
                minimizerList += myDir + " "
            elif mode == "density":
                densityList += myDir + " "
            elif mode == "genetic":
                geneticList += myDir + " "
            elif mode == "onerun":
                onerunList += myDir + " "
            elif mode == "diffusion":
                diffusionList += myDir + " "
                distributionList += myDir + " "
            else:
                print "Unknown mode"
    if plotAll or "minimizer" in argv:
        subprocess.call("python minimizerplot.py " + minimizerList, shell=True)
    if plotAll or "density" in argv:
        subprocess.call("python densityplot.py " + densityList, shell=True)
    if plotAll or "genetic" in argv:
        subprocess.call("python geneticplot.py " + geneticList, shell=True)
    if plotAll or "blocking" in argv:
        subprocess.call("python blockingplot.py " + onerunList, shell=True)
    if plotAll or "diffusion" in argv:
        subprocess.call("python diffusionplot.py " + diffusionList, shell=True)
    if plotAll or "distribution" in argv:
        subprocess.call("python distributionplot.py " + distributionList, shell=True)
print "Automator finished."
