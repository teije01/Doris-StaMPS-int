#!/usr/bin/env python
import numpy as np
import sys

arg = sys.argv;

#print arg

#for ar in arg[1:len(arg)]:
	#print ar
	
if len(arg) < 9:
	print "Needs at least eight number of arguments"
	print "Usage:"
	print "python  CropData.py 'orig_data_loc'  ori_lines  ori_pixels  'crop_data_loc'  line_0  line_1  pix_0  pix_1 [optional] dt "
	sys.exit("Wrong number of input arguments")

def CropData( orig_data_loc, ori_lines, ori_pixels, crop_data_loc, line_0, line_1, pix_0, pix_1, dt = np.dtype(np.complex64)):
	no_lines = (line_1 - line_0) + 1
	no_pixels = (pix_1 - pix_0) + 1
	crop = np.memmap(crop_data_loc, dtype=dt, mode='w+', shape=(no_lines, no_pixels))
	full_image = np.memmap(orig_data_loc, dtype=dt, mode='r', shape=(ori_lines, ori_pixels)) #order C
	crop[:no_lines,:no_pixels] = full_image[(line_0-1):line_1,(pix_0-1):pix_1]
	
def datatypefun( dattyp ):
	if dattyp in ['complex64', 'cr4', '8b']:
		datatype = 'complex64'
	elif dattyp in ['complex32', 'float32', '4b']:
		datatype = 'float32'
	else:
		print "Datatype unknown, please enter correct datatype"
		sys.exit("Create new functionality")
	dt = np.dtype( datatype )
	return dt
	
	
orig_data   = arg[1]
orig_lines  = int(float(arg[2]))
orig_pixels = int(float(arg[3]))
crop_data = arg[4]
l0 = int(float(arg[5]))
l1 = int(float(arg[6]))
p0 = int(float(arg[7]))
p1 = int(float(arg[8]))	

if len(arg) == 9:
	print "Cropping file: ", orig_data, "to file: ",crop_data
	CropData( orig_data, orig_lines, orig_pixels, crop_data, l0, l1, p0, p1 )
elif len(arg) == 10:
	dt = datatypefun( arg[9] )
	print "Cropping file: ", orig_data, "to file: ",crop_data
	CropData( orig_data, orig_lines, orig_pixels, crop_data, l0, l1, p0, p1, dt )
else:
	print "Not functional yet: other data type than Complex64 or Complex32"
	sys.exit("Create new functionality")


