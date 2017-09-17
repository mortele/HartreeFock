from overlapIntegralTests import *


if __name__ == '__main__' :
	for i in range(10) :
		for j in range(10) :
			Gb = G(*primitives[j])
			product = G(*primitives[i]) * (diff(Gb,x,x) + diff(Gb,y,y) + diff(Gb,z,z))
			integral = N(integrate(integrate(integrate(product, 
						(x, -oo, oo)), 
						(y, -oo, oo)), 
						(z, -oo, oo)))
			print("%d,%d " % (i,j), -0.5*integral)		