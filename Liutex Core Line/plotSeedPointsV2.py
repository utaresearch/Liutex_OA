#-------------------------------------------------------------------------------
# Plots Liutex Core Lines in Tecplot360 using seed points generated 
# by Get_Seeds_Liu_Core_Line.f95
#
# Author: Oscar Alvarez (oscar.alvarez@mavs.uta.edu)
# Summer 2022
#-------------------------------------------------------------------------------

import logging
import os
import tecplot as tp
from tecplot.constant import *
import numpy as np
import sys


logging.basicConfig(level=logging.DEBUG)

print("Program Starting")

if '-c' in sys.argv:
	print("Opening tp")
	tp.session.connect()

# Read input file
with open('input.txt', 'r') as f:
	f_start = int(f.readline().split()[0])
	f_end = int(f.readline().split()[0])
	f_skip = int(f.readline().split()[0])
	f.readline()
	f.readline()
	f.readline()
	f.readline()
	data_file_prefix = f.readline().split()[0]

# Find home directory (folder) where files are located
path = os.getcwd()

for i in range(f_start, f_end + f_skip, f_skip):

	tp.new_layout()
	
	grid_filename = data_file_prefix + '_' + str(i) + '.xyz'
	func_filename = 'seeds_' + data_file_prefix + '_' + str(i) + '.fun'
	seed_filename = 'seeds_' + data_file_prefix + '_' + str(i) + '.txt'
	
	# Read seed points file
	s_points = np.array([])
	
	with open(seed_filename, 'r') as f:
		for line in f:
			s_points = np.append(s_points, np.array(line.split()))

	print(seed_filename + " Read")

	l_s = len(s_points)
	
	print("Found " + str(int(l_s/3)) + " seeds")
	
	if(l_s == 0):
		print("!! Seed file is empty! NO SEED POINTS FOUND.")
		exit()

	# Load grid and function files
	grid_file = os.path.join(path, grid_filename)
	func_file = os.path.join(path, func_filename)
	seed_file = os.path.join(path, seed_filename)

	dataset = tp.data.load_plot3d(grid_filenames=grid_file, function_filenames=func_file)
	
	print(data_file_prefix + '_' + str(i) + " Read.")
	
	## Configuring Settings

	# Setup frame as Cartesinan 3D
	frame = tp.active_frame()
	frame.plot_type = tp.constant.PlotType.Cartesian3D
	
	plot = frame.plot()

	## Set up Variable names
	# Rename variables
	frame.dataset.variable('F1V8').name = "Liutex x"
	frame.dataset.variable('F1V9').name = "Liutex y"
	frame.dataset.variable('F1V10').name = "Liutex z"
	frame.dataset.variable('F1V11').name = "Liutex Magnitude"
	# frame.dataset.variable('F1V18').name = "Seed Point"

	# Set up stream trace plotting variables
	plot.vector.u_variable = dataset.variable('Liutex x')   # Liutex_x
	plot.vector.v_variable = dataset.variable('Liutex y')   # Liutex_y
	plot.vector.w_variable = dataset.variable('Liutex z')	# Liutex_z

	# Streamtrace settings
	streamtraces = plot.streamtraces
	streamtraces.color = Color.MultiColor2	# Multicolor to show Liutex Core Line
	streamtraces.show_arrows = True
	streamtraces.line_thickness = 0.2

	# Iso-Surface settings
	# plot.show_isosurfaces = True

	# iso = plot.isosurface(0)
	# iso.show = True
	# iso.definition_contour_group = plot.contour(0)
	# iso.isosurface_selection = IsoSurfaceSelection.OneSpecificValue
	# iso.isosurface_values = 1

	# plot.contour(1).variable = dataset.variable('Seed Point')
	# iso.contour.show = True
	# iso.contour.contour_type = ContourType.PrimaryValue
	# iso.contour.flood_contour_group = plot.contour(1)

	# Setting Liutex magnitude as the streamtrace color variable
	plot.contour(0).variable = dataset.variable('Liutex Magnitude')	# Liutex Magnitude for Color
	
	plot.contour(0).levels.reset_levels(np.linspace(0,.3,10))
	plot.show_contour = True
	plot.fieldmap(0).contour.show = True
	
	legend = plot.contour(0).legend
	legend.position = (90, 40) # Frame percentages

	#streamtraces.arrowhead_size = 1

	## Top View Settings
	# plot.view.width = 130	
	# plot.view.alpha = 0
	# plot.view.theta = 90
	# plot.view.psi   = 0
	# plot.view.position = (435.659, 14.1247, 844.763)

	## Plotting streamtraces and Positioning data for the picture

	print("Plotting Liutex Core Lines...")
	
	for j in range(0, l_s, 3):			
		# Plotting seed point
		seed_point = (s_points[j], s_points[j+1], s_points[j+2])
		stream_type = Streamtrace.VolumeLine
		streamtraces.add(seed_point, stream_type)


	# Save .png image
	img_filename = data_file_prefix + '_' + str(i) + '_image.png'
	img_file = os.path.join(path, img_filename)
	tp.export.save_png(img_file, 4000, supersample=2)
	print(img_filename + " Image Saved")

	# Save Tecplot Package
	layout_filename = data_file_prefix + '_' + str(i) + '_layout.lpk'
	tp.save_layout(layout_filename)
	print(layout_filename + " Tecplot Layout Saved")

