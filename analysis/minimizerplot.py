# -*- coding: utf-8 -*-
try:
    from enthought.mayavi.mlab import *
except ImportError:
    from mayavi.mlab import *

from numpy import *
from matplotlib import rc
from sys import argv
options.offscreen = True
first = True

for datapath in argv:
    if first:
	first = False
	continue
    try:
        param0 = loadtxt(datapath + "/parameters0.dat")
        param1 = loadtxt(datapath + "/parameters1.dat")
        data = loadtxt(datapath + "/energies.dat")
    except IOError as e:
        print "Could not find data files for " + myDir + ". Skipping..."
        continue
    print "Plotting " + datapath

    fig = figure(size=(960,960))
    fig.scene.background = (1.0, 1.0, 1.0)
    fig.scene.foreground = (0.0, 0.0, 0.0)

    surface = surf(param0, param1, data, warp_scale='auto')
    surface.actor.property.edge_visibility = True
    surface.actor.property.lighting = False
    
    myoutline = outline()

    myaxes = axes(xlabel='a', ylabel='b', zlabel='E')
    myaxes.title_text_property.italic = False
    myaxes.title_text_property.bold = False
    myaxes.title_text_property.font_family = 'times'
    myaxes.label_text_property.italic = False
    myaxes.label_text_property.bold = False
    myaxes.label_text_property.font_family = 'times'
    myaxes.axes.use_ranges = True
    myaxes.axes.corner_offset = 0.02
    
    # perspective view
    fig.scene.camera.position = [21.583470159040861, -16.378051020433883, 10.231416057474506]
    fig.scene.camera.focal_point = [5.3897238323904126, 5.434549130831372, -1.7158449313567394]
    fig.scene.render()

    savefig(datapath + "/minimizer.png")
    
    # alpha view
    fig.scene.x_plus_view()
    savefig(datapath + "/minimizer-alpha.png")
    
    # alpha view
    fig.scene.x_plus_view()
    fig.scene.parallel_projection = True
    myaxes.axes.y_axis_visibility = False
    myaxes.axes.z_axis_visibility = True
    savefig(datapath + "/minimizer-alpha.png")
    
    # beta view
    fig.scene.y_plus_view()
    fig.scene.camera.set_roll(180)
    fig.scene.parallel_projection = True
    myaxes.axes.x_axis_visibility = True
    myaxes.axes.y_axis_visibility = True
    myaxes.axes.z_axis_visibility = False
    savefig(datapath + "/minimizer-beta.png")