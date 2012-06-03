# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = '16'
import subprocess
import scipy.ndimage as ndimage
from scipy.interpolate import interp1d
from numpy import *
from matplotlib import rc
from sys import argv
import ConfigParser
from numpy.linalg.linalg import norm

first = True

for datapath in argv:
    try:
        from enthought.mayavi.mlab import *
    except ImportError:
        from mayavi.mlab import *
    #options.offscreen = True
    
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
    
    if len(data.shape) == 1:
	print "Fixing radial plot!"
	f2 = interp1d(param0, data, kind='cubic')
	param = concatenate((-param0[::-1], param0))
	mygrid = meshgrid(param,param)
	radius = sqrt(mygrid[0]**2 + mygrid[1]**2)
	radius = radius + (param0.max() - radius) * ( radius > param0.max())
	param0old = param0
	param1old = param1
	dataold = data
	param0 = mygrid[1]
	param1 = mygrid[0]
	data = f2(radius)
    
    fig = figure(size=(1280,960))

    fig.scene.background = (1.0, 1.0, 1.0)
    fig.scene.foreground = (0.0, 0.0, 0.0)

    theNorm = norm(data)
    if theNorm == inf:
	theNorm = data.max()
    if theNorm == inf:
	print "Found infinite data value. Skipping..."
        continue
    normalizedData = data / theNorm
    
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
            fig.scene.camera.position = [-69.07717592751338, -43.634636030649105, 42.06008568491923]
            fig.scene.camera.focal_point = [5.6484023160060648, 7.1971704691513159, -26.413896160865278]     
        else:
            fig.scene.camera.position = [-52.180466454173839, -19.827149342881295, 37.613833505235903]
            fig.scene.camera.focal_point = [11.309528075254283, 18.824854213069088, -7.9780393346049383]
    elif int(nParticles) == 6:
        aavava=0
    #    fig.scene.camera.position = [-83.367976780452835, -33.31347241689658, 60.149674418839126]
    #    fig.scene.camera.focal_point = [9.5877242100828539, 23.27692598937039, -6.6013866059718787]
    fig.scene.render()
    savefig(datapath + "/density.png")
    subprocess.call("convert " + datapath + "/density.png -trim " + datapath + "/density-trim.png", shell=True)
    
    # plot 2D
    from pylab import *
    figure()
    imshow(normalizedDataSmooth, extent=(param0.min(), param0.max(), param1.min(), param1.max()), origin="lower", interpolation='bilinear')
    cb = colorbar(format='$%g$')
    
    contour(param0, param1, normalizedDataSmooth, colors="k")
    ax = axes()
    formatter = mpl.ticker.FormatStrFormatter('$%g$')
    ax.xaxis.set_major_formatter(formatter) 
    ax.yaxis.set_major_formatter(formatter)
    savefig(datapath + "/density2d.pdf")
