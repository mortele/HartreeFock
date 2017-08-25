import numpy as np
from sympy.abc import *
from sympy import *


def G(i,j,k,a,Ax,Ay,Az,r):
	A = np.array([Ax,Ay,Az])
	rA = r - A
	Ax = rA[0]
	Ay = rA[1]
	Az = rA[2]
	rA2 = np.dot(rA,rA)
	return Ax**i * Ay**j * Az**k * exp(-a * rA2)


r = np.array([x, y, z])

primitive0 = [0,0,0,	1.0,	0.0, -1.5,  1.2,	r]
primitive1 = [1,0,0,	1.1,	0.1, -1.2,  2.5,	r]
primitive2 = [0,1,0,	0.3,	0.2, -1.0, -0.3,	r]
primitive3 = [0,0,1,	0.9,	0.3, -0.8,  0.1,	r]
primitive4 = [0,0,2,	2.2,	0.4, -0.6, -3.1,	r]
primitive5 = [0,2,0,	2.4,	0.5, -0.4,  3.8,	r]
primitive6 = [2,0,0,	3.1,	0.6, -0.2,  1.3,	r]
primitive7 = [1,1,0,	3.7,	0.7,  0.0,  2.4,	r]
primitive8 = [0,1,1,	1.3,	0.8,  0.2,  5.3,	r]
primitive9 = [1,0,1,	4.3,	0.9,  0.4,  1.2,	r]

primitives = [	primitive0,
				primitive1,
				primitive2,
				primitive3,
				primitive4,
				primitive5,
				primitive6,
				primitive7,
				primitive8,
				primitive9 ]

if __name__ == '__main__' :
	for i in range(10) :
		for j in range(10) :
			product = G(*primitives[i])*G(*primitives[j])
			integral = N(integrate(integrate(integrate(product, 
						(x, -oo, oo)), 
						(y, -oo, oo)), 
						(z, -oo, oo)))
			print("%d,%d " % (i,j), integral)
