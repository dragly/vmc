# -*- coding: utf-8 -*-
try:
    from enthought.mayavi.mlab import *
except ImportError:
    from mayavi.mlab import *

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
    surface = surf(param0, param1, normalizedData, warp_scale='auto')
    surface.actor.property.edge_visibility = True
    surface.actor.property.lighting = True

    myaxes = axes(xlabel='x', ylabel='y', zlabel='P(x,y)')
    myaxes.title_text_property.italic = False
    myaxes.title_text_property.bold = False
    myaxes.title_text_property.font_family = 'times'
    myaxes.label_text_property.italic = False
    myaxes.label_text_property.bold = False
    myaxes.label_text_property.font_family = 'times'
    myaxes.axes.use_ranges = True
    myaxes.axes.corner_offset = 0.1
    
    if nParticles == 2:
        fig.scene.camera.position = [73.692558727450901, 74.463726085840136, 33.595723834665613]
        fig.scene.camera.focal_point = [16.746857621650065, 17.518024980039424, -4.3680769025348729]
    fig.scene.render()

    savefig(datapath + "/density.png")