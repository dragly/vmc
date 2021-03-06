# -*- coding: utf-8 -*-
from numpy import *
from matplotlib import rc
from sys import argv
import subprocess
import matplotlib as mpl
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = '16'
first = True

for datapath in argv:
    if first:
	first = False
	continue
    try:
	from enthought.mayavi.mlab import *
    except ImportError:
	from mayavi.mlab import *
    options.offscreen = True
    print "Plotting " + datapath
    try:
        param0 = loadtxt(datapath + "/parameters0.dat")
        param1 = loadtxt(datapath + "/parameters1.dat")
        data = loadtxt(datapath + "/energies.dat")
    except IOError as e:
        print "Could not find data files for " + myDir + ". Skipping..."
        continue

    fig = figure(size=(960,960))
    fig.scene.background = (1.0, 1.0, 1.0)
    fig.scene.foreground = (0.0, 0.0, 0.0)

    surface = surf(param0, param1, data, warp_scale='auto')
    surface.actor.property.edge_visibility = True
    surface.actor.property.lighting = False
    
    myoutline = outline()

    myaxes = axes(xlabel=r'alpha', ylabel=r'beta', zlabel='E')
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
    
    # perspective view
    #fig.scene.camera.position = [25, -20, 12]
    #fig.scene.camera.focal_point = [5.3897238323904126, 5.434549130831372, -1.7158449313567394]
    #fig.scene.render()

    savefig(datapath + "/minimizer.png")
    subprocess.call("convert " + datapath + "/minimizer.png -trim " + datapath + "/minimizer-trim.png", shell=True)
    
    # alpha view
    #fig.scene.x_plus_view()
    #savefig(datapath + "/minimizer-alpha.png")
    
    # alpha view
    fig.scene.x_minus_view()
    fig.scene.parallel_projection = True
    myaxes.axes.y_axis_visibility = False
    myaxes.axes.z_axis_visibility = True
    savefig(datapath + "/minimizer-alpha.png")
    subprocess.call("convert " + datapath + "/minimizer-beta.png -trim " + datapath + "/minimizer-alpha-trim.png", shell=True)
    
    # beta view
    fig.scene.y_minus_view()
    fig.scene.camera.set_roll(0)
    fig.scene.parallel_projection = True
    myaxes.axes.x_axis_visibility = True
    myaxes.axes.y_axis_visibility = True
    myaxes.axes.z_axis_visibility = False
    savefig(datapath + "/minimizer-beta.png")
    subprocess.call("convert " + datapath + "/minimizer-alpha.png -trim " + datapath + "/minimizer-beta-trim.png", shell=True)
    
    # plot 2D
    from pylab import *
    figure()
    ax = axes([0.15, 0.12, 0.8, 0.83]) 
    imshow(transpose(data), extent=(param0.min(), param0.max(), param1.min(), param1.max()), origin="lower", interpolation='bilinear')
    cb = colorbar(format='$%g$')
    title(r"Energy")
    contour(param0, param1, data, 20, colors="k")
    formatter = mpl.ticker.FormatStrFormatter('$%g$')
    ax.xaxis.set_major_formatter(formatter) 
    ax.yaxis.set_major_formatter(formatter)
    xlabel(r"$\alpha$")
    ylabel(r"$\beta$")
    savefig(datapath + "/minimizer2d.pdf")
