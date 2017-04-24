#!/usr/bin/env python
import numpy as np
import sys
import numpy.polynomial.polynomial as poly
from scipy import interpolate
from scipy.optimize import fmin

arg = sys.argv;

#print arg

def readorb(filename):
    with open(filename, "r") as ins:
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

def orbit_int(orbits, method = 'poly', degree = 4):
    t = orbits[:,0]
    x = orbits[:,1]
    y = orbits[:,2]
    z = orbits[:,3]
    if method == 'poly':
        print("Using: poly -", degree)
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

def findtmin(xyz, tstart, coef, method = 'poly'):
    L = len(xyz)
    tmin = np.zeros((L,1))
    for idx in range(0,L):
        tmin[idx] = fmin(SR_dist, tstart, args=(coef, xyz[idx:idx+1,:], method), 
                    xtol=1e-10, ftol=1e-10, maxiter=250, maxfun=500, disp=0)
    return tmin

def Baseline(xyz_m, xyz_s):
    L = len(xyz_m)
    B = np.zeros((L,1))
    for idx in range(0,L):
        B[idx,0] = np.sqrt( (xyz_m[idx,0] - xyz_s[idx,0])**2 + 
                            (xyz_m[idx,1] - xyz_s[idx,1])**2 +
                            (xyz_m[idx,2] - xyz_s[idx,2])**2 )
    return B

def ParBaseline(xyz, xyz_m, xyz_s, step = 100):
    L  = len(xyz)
    Bp = np.zeros((L,1))
    for idx in range(0,L):
        idxs = int(np.floor(idx/step))
        dmp = np.sqrt( (xyz_m[idxs,0] - xyz[idx,0])**2 + 
                       (xyz_m[idxs,1] - xyz[idx,1])**2 +
                       (xyz_m[idxs,2] - xyz[idx,2])**2 )
        dsp = np.sqrt( (xyz_s[idxs,0] - xyz[idx,0])**2 + 
                       (xyz_s[idxs,1] - xyz[idx,1])**2 +
                       (xyz_s[idxs,2] - xyz[idx,2])**2 )
        Bp[idx,0] = dmp - dsp
    return Bp

def PerBaseline(Bpar, B, sign, step = 100):
    L  = len(Bpar)
    Bp = np.zeros((L,1))
    for idx in range(0,len(Bpar)):
        idxs = int(np.floor(idx/step))
        Bp[idx,0] = sign * np.sqrt( np.abs(B[idxs,0]**2 - Bpar[idx,0]**2) )
    return Bp

def detsign(xyz, xyz_m, xyz_s):
    P = xyz[50,:]
    M = xyz_m[0,:]
    S = xyz_s[0,:]
    r1 = M-P
    r2 = S-P
    a2 = np.sum(np.square(P));    a   = np.sqrt(a2)
    b2_1 = np.sum(np.square(r1)); b_1 = np.sqrt(b2_1)
    b2_2 = np.sum(np.square(r2)); b_2 = np.sqrt(b2_2)
    c2_1 = np.sum(np.square(M));
    c2_2 = np.sum(np.square(S)); 
    gam_1 = np.arccos((a2+b2_1-c2_1)/(2*a*b_1)) * 180 / np.pi
    gam_2 = np.arccos((a2+b2_2-c2_2)/(2*a*b_2)) * 180 / np.pi
    #print(gam_1,gam_2)
    if gam_1 < gam_2:
        sign = 1
    elif gam_1 >= gam_2:
        sign = -1
    else:
        print("error")
    return sign
	
def Bhva(Bper, Bpar, theta):
	theta_r = theta * np.pi / 180
	Bh = (
	np.multiply(Bper[:,0], np.cos(theta_r)) + 
	np.multiply(Bpar[:,0], np.sin(theta_r))  )
	Bv = (
	np.multiply(Bper[:,0], np.sin(theta_r)) - 
	np.multiply(Bpar[:,0], np.cos(theta_r))  )
	alpha = (np.arctan2(Bv,Bh))
	alpha_d = alpha * 180 / np.pi
	return Bh,Bv,alpha_d
	
if len(arg) < 2:
	print "Needs at least one argument"
	print "Usage:"
	print "python  BperpDate.py   'Filename'  method   degree"
	print "If method choose: 'poly' or 'spline' with an uneven degree for splines"
	sys.exit("Wrong number of input arguments")

FN = arg[1]

if len(arg) == 2:
	met = 'poly'
	deg = 4
elif len(arg) == 3:
	met = arg[2]
	deg = 4
elif len(arg) == 4:
	met = arg[2]
	deg = float(arg[3])
else:
	print "Not functional yet: more than 3 input arguments"
	sys.exit("Create new functionality")

print ""
print "-- Starting python script: BperpDate.py --"
print "Reading master and slave orbits"
morbits = readorb('master.res')
sorbits = readorb('slave.res' )
print "Importing master azimuth time and data xyz points"
theta   = np.loadtxt('../look_angle.1.in')
taz_m   = np.loadtxt('../tmin.txt')
xyz     = np.loadtxt('../xyz.txt')
print "Interpollation of master and slave orbits"
mcoef = orbit_int(morbits, method = met, degree = deg)
scoef = orbit_int(sorbits, method = met, degree = deg)
xyz_m = xyz_t(mcoef, taz_m, method = met)
xyz_m = np.transpose(xyz_m)
print "Calculating slave azimuth time"
taz_s = findtmin(xyz_m, np.mean(taz_m), scoef, method = met)
taz_s = np.transpose(taz_s)
xyz_s = xyz_t(scoef, taz_s, method = met)
xyz_s = np.transpose(xyz_s)
print "Calculating Baseline parameters:"
B    = Baseline(xyz_m, xyz_s)
Bpar = ParBaseline(xyz, xyz_m, xyz_s)
sign = detsign(xyz, xyz_m, xyz_s)
Bper = PerBaseline(Bpar, B, sign)
Bh, Bv, a = Bhva(Bper, Bpar, theta)
Bpars = np.array([np.mean(B), np.mean(Bpar), np.mean(Bper), 
				   np.mean(Bh), np.mean(Bv), np.mean(a)])
print "B    : ",Bpars[0]
print "Bpar : ",Bpars[1]
print "Bper : ",Bpars[2]
print "Bh   : ",Bpars[3]
print "Bv   : ",Bpars[4]
print "alpha: ",Bpars[5]
print "Saving Bperp parameter file in", FN
np.savetxt( FN, Bper, fmt='%6.4f' )
np.savetxt( "Baseline.pars", Bpars, fmt='%6.6f')
print "-- Python script finished --"
print ""

