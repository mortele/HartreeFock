using Cubature


const p =	[0 0 0 1.0 0.0 -1.5  1.2;
		   	 1 0 0 1.1 0.1 -1.2  2.5;
		   	 0 1 0 0.3 0.2 -1.0 -0.3;
		   	 0 0 1 0.9 0.3 -0.8  0.1;
		   	 0 0 2 2.2 0.4 -0.6 -3.1;
		   	 0 2 0 2.4 0.5 -0.4  3.8;
		   	 2 0 0 3.1 0.6 -0.2  1.3;
		   	 1 1 0 3.7 0.7  0.0  2.4;
		   	 0 1 1 1.3 0.8  0.2  5.3;
		   	 1 0 1 4.3 0.9  0.4  1.2]


function G( i ::Int, 	j ::Int, 	k ::Int, 	a::Number, 
			Ax::Number, Ay::Number, Az::Number, 
			x ::Number, y ::Number, z ::Number)
	rAx = x-Ax
	rAy = y-Ay
	rAz = z-Az
	rA2 = rAx^2 + rAy^2 + rAz^2
	return rAx^i * rAy^j * rAz^k * exp(-a * rA2)
end

function Gp( i1 ::Int, 		j1 ::Int, 		k1 ::Int, 		a1::Number, 
			 Ax1::Number, 	Ay1::Number, 	Az1::Number, 
			 i2 ::Int, 		j2 ::Int, 		k2 ::Int, 		a2::Number, 
			 Ax2::Number, 	Ay2::Number, 	Az2::Number,
			 x  ::Number, 	y  ::Number, 	z  ::Number,
			 nx ::Number, 	ny ::Number, 	nz ::Number)
	rA1x = x-Ax1
	rA1y = y-Ay1
	rA1z = z-Az1
	rA12 = rA1x^2 + rA1y^2 + rA1z^2
	G1   = rA1x^i1 * rA1y^j1 * rA1z^k1 * exp(-a1 * rA12)

	rA2x = x-Ax2
	rA2y = y-Ay2
	rA2z = z-Az2
	rA22 = rA2x^2 + rA2y^2 + rA2z^2
	G2   = rA2x^i2 * rA2y^j2 * rA2z^k2 * exp(-a2 * rA22)

	rCx  = x-nx
	rCy  = y-ny
	rCz  = z-nz
	rC   = sqrt(rCx^2 + rCy^2 + rCz^2)

	return G1 * (1.0/rC) * G2
end

const limit = 30.0

for i_ = 1:10
	for j_ = 1:10
		iii   = i_-1
		jjj   = j_-1
		atomm = (359*iii+295*jjj-42*iii*jjj+120*iii*iii-38*iii*jjj*iii) % 9
        atomm = (atomm > 0 ? atomm : -atomm)
        print(atomm, " ")

		function integrand(xv::Array{Float64,1})
			x = xv[1]
			y = xv[2]
			z = xv[3]

			ii   = i_-1
			jj   = j_-1
			atom = (359*ii+295*jj-42*ii*jj+120*ii*ii-38*ii*jj*ii) % 9
            atom = (atom > 0 ? atom : -atom)
            nx   = p[atom+1,5]
            ny   = p[atom+1,6]
            nz   = p[atom+1,7]

			i1  = Int(p[i_,1])
			j1  = Int(p[i_,2])
			k1  = Int(p[i_,3])
			a1  = p[i_,4]
			Ax1 = p[i_,5]
			Ay1 = p[i_,6]
			Az1 = p[i_,7]

			i2  = Int(p[j_,1])
			j2  = Int(p[j_,2])
			k2  = Int(p[j_,3])
			a2  = p[j_,4]
			Ax2 = p[j_,5]
			Ay2 = p[j_,6]
			Az2 = p[j_,7]

			return Gp(	i1,j1,k1,a1,
						Ax1,Ay1,Az1, 
					 	i2,j2,k2,a2, 
					 	Ax2,Ay2,Az2,
					 	x,y,z,
					 	nx,ny,nz)
		end

		(val,err) = hcubature(integrand, 	[-limit -limit -limit], 
											[ limit  limit  limit],
											maxevals=Int(5e8))
        println(i_-1,",",j_-1,"  ",val, "  ", err)
	end
end