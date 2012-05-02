from enthought.mayavi.mlab import *
from numpy import *
from matplotlib import rc
param0 = loadtxt("../../vmc-release/parameters0.dat")
param1 = loadtxt("../../vmc-release/parameters1.dat")
data = loadtxt("../../vmc-release/energies.dat")

surface = surf(param0, param1, data, warp_scale='auto')
surface.actor.property.edge_visibility = True

axes = axes(xlabel=r'$\alpha$', ylabel=r'$\beta$', zlabel=r'$E$')
axes.title_text_property.italic = False
axes.title_text_property.bold = False
axes.title_text_property.font_family = 'times'
axes.label_text_property.italic = False
axes.label_text_property.bold = False
axes.label_text_property.font_family = 'times'
axes.axes.use_ranges = True
