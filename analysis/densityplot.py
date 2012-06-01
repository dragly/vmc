# -*- coding: utf-8 -*-
try:
    from enthought.mayavi.mlab import *
except ImportError:
    from mayavi.mlab import *
import subprocess
import scipy.ndimage as ndimage
from numpy import *
from matplotlib import rc
from sys import argv
import ConfigParser
options.offscreen = True
first = True

for datapath in argv:
    if first:
	first = False
	continue

    try:
        config = ConfigParser.ConfigParser()
        configFileName = datapath + "/config.ini"
        config.read(configFileName)
        nParticles = config.get("General", "nParticles")
        omega = config.get("General", "omega")

        param0 = loadtxt(datapath + "/density-grid-x.dat")
        param1 = loadtxt(datapath + "/density-grid-y.dat")
        data = loadtxt(datapath + "/density.dat")
    except IOError as e:
        print "Could not find data files for " + myDir + ". Skipping..."
        continue
    print "Plotting " + datapath
    
    fig = figure(size=(1280,960))

    fig.scene.background = (1.0, 1.0, 1.0)
    fig.scene.foreground = (0.0, 0.0, 0.0)

    from numpy.linalg.linalg import norm
    normalizedData = data / norm(data)
    
    normalizedDataSmooth = ndimage.gaussian_filter(normalizedData, sigma=1.0, order=0)
    
    surface = surf(param0, param1, normalizedDataSmooth, warp_scale='auto')
    surface.actor.property.edge_visibility = True
    surface.actor.property.lighting = True
    #outline()

    myaxes = axes(xlabel='x', ylabel='y', zlabel='P(x,y)')
    myaxes.title_text_property.italic = False
    myaxes.title_text_property.bold = False
    myaxes.title_text_property.font_family = 'times'
    myaxes.label_text_property.italic = False
    myaxes.label_text_property.bold = False
    myaxes.label_text_property.font_family = 'times'
    myaxes.axes.use_ranges = True
    myaxes.axes.corner_offset = 0.1
    myaxes.axes.number_of_labels = 4
    myaxes.axes.label_format = '%-#6.5g'
    
    if int(nParticles) == 2:
        if float(omega) < 0.20:
            fig.scene.camera.position = [-67.71549063812715, -40.106253119177815, 46.165393868709053]
            fig.scene.camera.focal_point = [7.0100876053923837, 10.725553380622589, -22.308587977075494]          
        else:
            fig.scene.camera.position = [-52.180466454173839, -19.827149342881295, 37.613833505235903]
            fig.scene.camera.focal_point = [11.309528075254283, 18.824854213069088, -7.9780393346049383]
    elif int(nParticles) == 6:
        fig.scene.camera.position = [-83.367976780452835, -33.31347241689658, 60.149674418839126]
        fig.scene.camera.focal_point = [9.5877242100828539, 23.27692598937039, -6.6013866059718787]
    fig.scene.render()
    savefig(datapath + "/density.png")
    subprocess.call("convert " + datapath + "/density.png -trim " + datapath + "/density-trim.png", shell=True)
    
    
