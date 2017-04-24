#!/usr/bin/env python
import numpy as np
import numpy.matlib
import sys
import numpy.polynomial.polynomial as poly
from scipy import interpolate
from scipy.optimize import fmin

arg = sys.argv;

#print arg
	
if len(arg) < 5:
	print "Needs at least four number of arguments"
	print "Usage:"
	print "python  LookAngle.py   line_0  line_1  pix_0  pix_1    method(opt)  degree(opt)"
	print "If method choose: 'poly' or 'spline' with an uneven degree for splines"
	sys.exit("Wrong number of input arguments")

def LPH( line_0, line_1, pix_0, pix_1 ):
	#pixels
	p = np.linspace(pix_0, pix_1, num=50)
	p = np.matlib.repmat(p,2,1)
	p = p.reshape(1,100,order='F')
	p = np.matlib.repmat(p,50,1)
	p = p.reshape(5000,1)
	#lines
	l = np.linspace(line_0, line_1 ,num=50)
	l = l.reshape(50,1)
	l = np.matlib.repmat(l,1,100)
	l = l.reshape(5000,1)
	#height
	h = np.array([0, 1000])
	h = np.matlib.repmat(h,1,2500)
	h = h.reshape(5000,1)
	#CAT
	lph = np.hstack((l,p,h))
	return lph
	
def LPH_corr( lph ):
	lph[:,0] = lph[:,0] - lph[0,0] + 1
	lph[:,1] = lph[:,1] - lph[0,1] + 1
	return lph
	
def RotMat(inp_vec, theta):
	TM = np.zeros((2,2))
	TM = [[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]]
	out_vec = np.dot(TM,inp_vec)
	return out_vec
	
def ext_pl( lph, dt='float32' ):
	plh = lph
	n_lines  = lph[-1,0]
	n_pixels = lph[-1,1]
	full_lam = np.memmap('lon.raw', dtype=np.dtype( dt ), mode='r', shape=(n_lines, n_pixels))
	full_phi = np.memmap('lat.raw', dtype=np.dtype( dt ), mode='r', shape=(n_lines, n_pixels))
	a = full_phi[-1,np.floor(n_pixels/2)] - full_phi[0,np.floor(n_pixels/2)]
	b = full_lam[-1,np.floor(n_pixels/2)] - full_lam[0,np.floor(n_pixels/2)]
	theta   = -(np.arctan2(a,b) - np.pi/2)
	#theta_d = -(np.arctan2(a,b) - np.pi/2) * (180 / np.pi)
	for pnt in range(0, len(lph)):
		l = lph[pnt,0]-1; lf = int(np.floor(l)); lc = int(np.ceil(l)) #index -1 because python
		p = lph[pnt,1]-1; pf = int(np.floor(p)); pc = int(np.ceil(p)) #index -1 because python
		if lc == lph[-1,0]-1:
			1+1
		elif lf == lc:
			lc = lc + 1
		if pc == lph[-1,1]-1:
			1+1
		elif pf == pc:
			pc = pc + 1
		dp = p - pf; dl = l - lf;
		#
		inp_vecb = np.array([dp,1-dl])
		[PH, b]    = RotMat(inp_vecb, theta)			#PH = placeholder
		[PH, b_sc] = RotMat(np.array([1,1]), theta)	#PH = placeholder
		B = b/b_sc;
		plh[pnt,0] = full_phi[lf,pc]*B + full_phi[lc,pf]*(1-B)
		#
		inp_veca = np.array([dl,dp])
		[PH, a]    = RotMat(inp_veca, theta)			#PH = placeholder
		[PH, a_sc] = RotMat(np.array([1,1]), theta)	#PH = placeholder
		A = a/a_sc;
		plh[pnt,1] = full_lam[lc,pc]*A + full_lam[lf,pf]*(1-A)
	return plh
	
def plh2xyz( plh, A=6378137.0, FL=0.0033528128981864433 ):
	#input parameters
	# * ----------------
	# * A                semi-major axis of ellipsoid [units are of distance]
	# * FL               flattening of ellipsoid [unitless]
	# * plh[]            ellipsoidal coordinates of point, in geodetic latitude,
	# *                  longitude east of Greenwich, and height [units for
	# *                  latitude=plh[0] and longitude=plh[1] are in degrees;
	# *                  height=plh[2] are distance and will be the same as
	# *                  those of the input parameters]
	flatfn = (2.0 - FL)*FL
	funsq  = (1.0 - FL)*(1.0 - FL)
	lat_rad = (np.pi/180) * plh[:,0]
	lon_rad = (np.pi/180) * plh[:,1]
	sin_lat = np.sin(lat_rad)
	#
	g1 = A / np.sqrt(1.0 - flatfn * sin_lat**2)
	g2 = g1*funsq + plh[:,2]
	g1 = g1 + plh[:,2]
	#
	xyz = np.zeros(np.shape(plh))
	xyz[:,0] = g1 * np.cos( lat_rad )
	xyz[:,1] = xyz[:,0] * np.sin( lon_rad )
	xyz[:,0] = xyz[:,0] * np.cos( lon_rad )
	xyz[:,2] = g2 * sin_lat
	return xyz
	
def readorb():
	with open("master.res", "r") as ins:
		array = []
		mode = "search"
		for line in ins:
			if mode == "search":
				if line.rstrip('\n') == '*_Start_precise_orbits:':
					mode = "read"
					skip = 3
				else:
					1+1
			elif mode == "read":
				if skip == 0:
					array.append(line.rstrip('\n'))
					if array[-1] == '******************************************************************* ':
						del array[-1]
						mode = "stop"
				elif skip > 0:
					skip = skip -1
			elif mode == "stop":
				1+1
		orbits = np.zeros((0,4))
		for line in array:
			arr = np.fromstring(line, dtype=float, sep=' ')
			orbits = np.vstack((orbits,arr[0:5]))
	return orbits
	
def orbit_int(orbits, method = 'poly', degree = 3):
    t = orbits[:,0]
    x = orbits[:,1]
    y = orbits[:,2]
    z = orbits[:,3]
    if method == 'poly':
        print "Using: poly -", degree
        x_coef = poly.polyfit(t, x, degree )
        y_coef = poly.polyfit(t, y, degree )
        z_coef = poly.polyfit(t, z, degree )
        coef = np.vstack((x_coef, y_coef, z_coef))
        return coef
    elif method == 'spline':
        print "Using: Spline -", degree
        if (float(degree)/2.0) == int(degree/2.0):
            print "It is not recommended to have an even degree for splines"
        smoothing = 3.0
        x_tck = interpolate.splrep(t, x, k=degree, s = smoothing)
        y_tck = interpolate.splrep(t, y, k=degree, s = smoothing)
        z_tck = interpolate.splrep(t, z, k=degree, s = smoothing)
        coef = {'x': x_tck, 'y': y_tck, 'z': z_tck } #dictionary
        return coef
    else:
        print "method unkown, change method"

def xyz_t(coef, t, method = 'poly'):
    if method == 'poly':
        x_fit = poly.polyval(t, coef[0,:])
        y_fit = poly.polyval(t, coef[1,:])
        z_fit = poly.polyval(t, coef[2,:])
    elif method == 'spline':
        x_tck = coef['x']
        y_tck = coef['y']
        z_tck = coef['z']
        x_fit = interpolate.splev(t, x_tck)
        y_fit = interpolate.splev(t, y_tck)
        z_fit = interpolate.splev(t, z_tck)
    else:
        print "method unkown, change method"
    fit = np.vstack((np.float64(x_fit), np.float64(y_fit), np.float64(z_fit)))
    return fit

def SR_dist(t, coef, xyz, method = 'poly'):
    #if np.shape(xyz[0:1,:]) == (1, 3):
    #    print('xyz in good format: xyz = ',xyz) 
    fit = xyz_t(coef, t, method)
    fit = np.double(fit)
    xyz = np.double(xyz)
    dist = np.sqrt( (fit[0,:] - xyz[0][0])**2 + (fit[1,:] - xyz[0][1])**2 +
                    (fit[2,:] - xyz[0][2])**2)
    return dist

def findtmin(xyz, tstart, coef, method = 'poly', step = 100):
    L = int(np.floor(len(xyz)/step))
    tmin = np.zeros((L,1))
    for idx in range(0,L):
        idxyz = idx*100
        tmin[idx] = fmin(SR_dist, tstart, args=(coef, xyz[idxyz:idxyz+1,:], method), 
                    xtol=1e-10, ftol=1e-10, maxiter=100, maxfun=200, disp=0)
    return tmin

def LookAngle(xyz, tmin, coef, method = 'poly', step = 100):
	theta = np.zeros((5000,1))
	incid = np.zeros((5000,1))
	for ind in range(0,5000):
		t_az = tmin[np.floor(ind/step)]
		p_orb = xyz_t(coef, t_az, method)
		p_orb = p_orb[0:3,0]  # point on orbit interpollation
		p_ear = xyz[ind]      # point on earth
		v_oe  = p_orb - p_ear # vector (orbit - earth)
		#
		a2 = np.sum(np.square(p_ear))
		a  = np.sqrt(a2)
		b2 = np.sum(np.square(v_oe))
		b  = np.sqrt(b2)
		c2 = np.sum(np.square(p_orb))
		c  = np.sqrt(c2)
		theta_rad = np.arccos((-a2+b2+c2)/(2*b*c))
		theta_deg = theta_rad * 180 / np.pi
		gamma_rad = np.arccos(( a2+b2-c2)/(2*a*b))
		gamma_deg = gamma_rad * 180 / np.pi
		incidence = 180 - gamma_deg
		theta[ind] = theta_deg
		incid[ind] = incidence
	return theta, incid
	
# Action section
	
l0 = int(float(arg[1]))
l1 = int(float(arg[2]))
p0 = int(float(arg[3]))
p1 = int(float(arg[4]))	

if len(arg) == 5:
	met = 'poly'
	deg = 4
elif len(arg) == 6:
	met = arg[5]
	deg = 5
elif len(arg) == 7:
	met = arg[5]
	deg = int(float(arg[6]))
else:
	print "Not functional yet: lower than 5 or higer than 7 input arguments"
	sys.exit("Create new functionality")

print ""
print "-- Starting python script: LookAngle.py --"
print "Create line pixel height raster for StaMPS: lph";
lph = LPH( l0, l1, p0, p1 )
print "Correct line pixel raster height";
lph = LPH_corr( lph )
print "Convert to lon lat height format";
plh = ext_pl( lph )
print "Convert to xyz format ";
xyz = plh2xyz( plh )
print "Read orbit from master.res";
orbits = readorb()
print "Interpollate orbit";
coef = orbit_int(orbits, method = met, degree = deg)
t = orbits[:,0]
tstart = np.mean((t[0],t[-1]))
print "Find satellite azimuth times"
tmin = findtmin(xyz, tstart, coef, method = met)
print "Calculate Look Angles and incidence angles for the master file"
theta, incid = LookAngle(xyz, tmin, coef, method = met)

FN1 = 'look_angle.1.in'
FN2 = 'incidence_angle.1.in'
FN3 = 'tmin.txt'
FN4 = 'xyz.txt'
print "Saving ",FN1
print "Saving ",FN2
print "Saving ",FN3
print "Saving ",FN4
np.savetxt( FN1, theta, fmt='%6.4f' )
np.savetxt( FN2, incid, fmt='%6.4f' )
np.savetxt( FN3, tmin, fmt='%6.10f' )
np.savetxt( FN4, xyz, fmt='%6.10f' )
print "-- Python script finished --"
print ""



