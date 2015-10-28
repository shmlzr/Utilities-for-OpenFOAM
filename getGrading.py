#!/usr/bin/env python
#	Calculation of grading factor in openFOAM blockMesh. 
#	Uses "Y+ wall distance estimation" according to http://www.cfd-online.com/Wiki/Y_plus_wall_distance_estimation 
#	
#
#	COMMAND LINE USE
#	python	 getGrading.py "U" "D" "yp" "L" "N" ("+/-") ("rho") ("nuKin")
#
#	
#	example:
#	python getGrading.py 1 0.044 1.0 0.022 100 -
#
#
#	INPUT
#	- U:	Bulk velocity
#	- D:	Pipe diameter
#	- yP:	Desired y+ value
#	- L: 	Geometric length in direction of grading in meter [m]
#	- N: 	Number of cells to be used
#
#     - +/-: Grading ratio. If grading is from small to large --> "-"
#	- rho:  Density of the fluid. Default: 1.205 [kg/m3]
#	- nuKin:Kinematic viscosity.  Default: 1e-6 [m2/s]

#	
#	OUTPUT
#	- Reynolds number according to the given boundary layer length
#	- smallest/biggest cell size
#	- grading factor to be used in blockMeshDict
#
#	
#     Martin Schmelzer,	
#	m.schmelzer@tudelft.nl
#


# python 
import sys
import subprocess
from math import *
import numpy as np

def float_round(num, places = 0, direction = floor):
    return direction(num * (10**places)) / float(10**places)

# declaration of input parameters
if (len(sys.argv) == 6):
	U	= float(sys.argv[1])	#[m/s] 	freestream velocity
	D	= float(sys.argv[2])	#[m]	boundary layer length
	yPlus	= float(sys.argv[3])	#[]	desired y+ 
	L 	= float(sys.argv[4])	#[m]	geometric length of domain along grading will be applied
	N 	= int(sys.argv[5])	#[]	number of cells along this domain
	gradRatio= "-"
	rho 	= 1.205 		#[kg/m3] density of air at T=20C
	nuKin = 1e-6 			#[m2/s]  kinematic viscosity
elif (len(sys.argv) == 7):
	U	= float(sys.argv[1])	#[m/s] 	freestream velocity
	D	= float(sys.argv[2])	#[m]	boundary layer length
	yPlus	= float(sys.argv[3])	#[]	desired y+ 
	L 	= float(sys.argv[4])	#[m]	geometric length of domain along grading will be applied
	N 	= int(sys.argv[5])	#[]	number of cells along this domain
	gradRatio= sys.argv[6]
	rho 	= 1.205			#[kg/m3]density of the fluid
	nuKin	= 1e-6			#[m2/s] kinematic viscosity of the fluid
elif (len(sys.argv) == 9):
	U	= float(sys.argv[1])	#[m/s] 	freestream velocity
	D	= float(sys.argv[2])	#[m]	boundary layer length
	yPlus	= float(sys.argv[3])	#[]	desired y+ 
	L 	= float(sys.argv[4])	#[m]	geometric length of domain along grading will be applied
	N 	= int(sys.argv[5])	#[]	number of cells along this domain
	gradRatio= sys.argv[6]
	rho 	= float(sys.argv[7])	#[kg/m3]density of the fluid
	nuKin	= float(sys.argv[8])	#[m2/s] kinematic viscosity of the fluid
elif (len(sys.argv) < 6 or len(sys.argv) > 9 or len(sys.argv) == 8):
	print "\n", "ERROR: You have to specify 5, 6 or 8 input parameters!\n "
	print "example:"
	print "python getGrading.py U D yp L N (+/-) (rho) (nuKin)\n", "python getGrading.py U D yp L N - rho nuKin\n"
	sys.exit()


# Calculation of smallest cell 
Re = U * D / nuKin			           #[] Reynolds number accordin to the pipe diameter
Cf = (2*log(Re)/log(10) - 0.65)**(-2.3)  #[] Schlichting skin friction for Re<10e9
tauW = Cf * 0.5 * rho * U**2		    #[kg/ms2] wall shear stress
uStar = sqrt(tauW / rho)		          #[m/s] friction velocity

cellSizeMin = (yPlus * nuKin) / (uStar) #[m] wall distance or smallest cell size


# Estimate constant grading ratio of neigbouring cells iteratively
itmax = 500
F = L / cellSizeMin
if (gradRatio=="+"):
	c = 1.05
	for i in range(0,itmax - 1):
		c = (F*(c-1.) + 1.)**(1./N)
	
	# Compute vector of cell sizes for whole domain
	dx =[]
	X=0.0
	for i in range(0,N):
		tmp = cellSizeMin*(c**i)
		dx.append(tmp)
		X=X+tmp

	# Compute grading factor for blockMeshDict
	Grading = round(c**(N-1.))
	
elif (gradRatio=="-"):
	c = 0.05
	for i in range(0,itmax - 1):
		c = (c - F*(c-1.) )**(1./(1-N))
	
	# Compute vector of cell sizes for whole domain
	dx =[]
	X=0.0
	for i in range(0,N):
		tmp = cellSizeMin*((1./c)**i)
		dx.append(tmp)
		X=X+tmp
	# Compute grading factor for blockMeshDict
	Grading = float_round(c**(N-1.), 3, round)


# Test the OpenFOAM implementation
Grading=c**(N-1.0)
r = Grading ** (1.0/(N-1.0))
if Grading>1.0:
    alpha = Grading
else:
    alpha = 1.0 - r**(-N) + r**(-1.0)
dxmin = L *(r-1.0)/(alpha *r -1.0)




# Time step according to CFL condition
CFL = np.array([0.3, 0.7, 1.0])
deltaT = CFL * cellSizeMin / U

	

print "\n"	
#print "Vector of cell size distribution	[m]: \n", dx, "\n"
#print "Geometric length [m]: \n", X , "\n"
print "Reynolds Number: \n", Re , "\n"
print "Desired y+ value: \n", yPlus , "\n"
print "Smallest cell size [m]: \n", cellSizeMin , "\n"
#print "dxmin nach Of \n", dxmin , "\n"
print "Biggest cell size [m]: \n", dx[N-1] , "\n"
print "Grading factor to be used in blockMeshDict: \n", Grading , "\n"

print "Time step for time-resolved simulations:"
for i in range(0,len(CFL)):
	print "According to CFL:", CFL[i], "--> dt =",deltaT[i], "[s]"

print "\n"

