# -*- coding: utf-8 -*-

#servers = "boblekammer ccd diode forsterker fotodiode frekvensteller geigerteller hodoskop hologram interferometer kalorimeter kondensator kretskort laser likning mikroprosessor mosfet motstand multimeter scintillator signal signalgenerator spole synkrotron termometer transduser transformator hartreefock gran"
servers = "boblekammer ccd diode forsterker fotodiode frekvensteller geigerteller hologram interferometer kalorimeter kondensator kretskort laser likning mikroprosessor mosfet motstand multimeter scintillator signal signalgenerator spole synkrotron termometer transduser transformator hartreefock gran"
runFolder = "runs"
configFileName = "config.ini"
prevrunConfigFileName = "config.ini.prevrun"
lockFileName = "lock"
runCommand = "./vmc"
runCommandMPI = "'mpirun -n 4 ./vmc'"