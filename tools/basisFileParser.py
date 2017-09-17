import sys 

names = {	"H"   : [1, "Hydrogen"],
			"He"  : [2, "Helium"],	
			"Li"  : [3, "Lithium"],	
			"Be"  : [4, "Beryllium"],	
			"B"   : [5, "Boron"],	
			"C"   : [6, "Carbon"],	
			"N"   : [7, "Nitrogen"],	
			"O"   : [8, "Oxygen"],	
			"F"   : [9, "Fluorine"],	
			"Ne"  : [10, "Neon"],	
			"Na"  : [11, "Sodium"],	
			"Mg"  : [12, "Magnesium"],	
			"Al"  : [13, "Aluminum"],	
			"Si"  : [14, "Silicon"],	
			"P"   : [15, "Phosphorus"],	
			"S"   : [16, "Sulfur"],
			"Cl"  : [17, "Chlorine"],
			"Ar"  : [18, "Argon"],
			"K"   : [19, "Potassium"],
			"Ca"  : [20, "Calcium"],
			"Sc"  : [21, "Scandium"],
			"Ti"  : [22, "Titanium"],
			"V"   : [23, "Vanadium"],
			"Cr"  : [24, "Chromium"],
			"Mn"  : [25, "Manganese"],
			"Fe"  : [26, "Iron"],
			"Co"  : [27, "Cobalt"],
			"Ni"  : [28, "Nickel"],
			"Cu"  : [29, "Copper"],
			"Zn"  : [30, "Zinc"],
			"Ga"  : [31, "Gallium"],
			"Ge"  : [32, "Germanium"],
			"As"  : [33, "Arsenic"],
			"Se"  : [34, "Selenium"],
			"Br"  : [35, "Bromine"],
			"Kr"  : [36, "Krypton"],
			"Rb"  : [37, "Rubidium"],
			"Sr"  : [38, "Strontium"],
			"Y"   : [39, "Yttrium"],
			"Zr"  : [40, "Zirconium"],
			"Nb"  : [41, "Niobium"],
			"Mo"  : [42, "Molybdenum"],
			"Tc"  : [43, "Technetium"],
			"Ru"  : [44, "Ruthenium"],
			"Rh"  : [45, "Rhodium"],
			"Pd"  : [46, "Palladium"],
			"Ag"  : [47, "Silver"],
			"Cd"  : [48, "Cadmium"],
			"In"  : [49, "Indium"],
			"Sn"  : [50, "Tin"],
			"Sb"  : [51, "Antimony"],
			"Te"  : [52, "Tellurium"],
			"I"   : [53, "Iodine"],
			"Xe"  : [54, "Xenon"],
			"Cs"  : [55, "Cesium"],
			"Ba"  : [56, "Barium"],
			"La"  : [57, "Lanthanum"],
			"Ce"  : [58, "Cerium"],
			"Pr"  : [59, "Praseodymium"],
			"Nd"  : [60, "Neodymium"],
			"Pm"  : [61, "Promethium"],
			"Sm"  : [62, "Samarium"],
			"Eu"  : [63, "Europium"],
			"Gd"  : [64, "Gadolinium"],
			"Tb"  : [65, "Terbium"],
			"Dy"  : [66, "Dysprosium"],
			"Ho"  : [67, "Holmium"],
			"Er"  : [68, "Erbium"],
			"Tm"  : [69, "Thulium"],
			"Yb"  : [70, "Ytterbium"],
			"Lu"  : [71, "Lutetium"],
			"Hf"  : [72, "Hafnium"],
			"Ta"  : [73, "Tantalum"],
			"W"   : [74, "Tungsten"],
			"Re"  : [75, "Rhenium"],
			"Os"  : [76, "Osmium"],
			"Ir"  : [77, "Iridium"],
			"Pt"  : [78, "Platinum"],
			"Au"  : [79, "Gold"],
			"Hg"  : [80, "Mercury"],
			"Tl"  : [81, "Thallium"],
			"Pb"  : [82, "Lead"],
			"Bi"  : [83, "Bismuth"],
			"Po"  : [84, "Polonium"],
			"At"  : [85, "Astatine"],
			"Rn"  : [86, "Radon"],
			"Fr"  : [87, "Francium"],
			"Ra"  : [88, "Radium"],
			"Ac"  : [89, "Actinium"],
			"Th"  : [90, "Thorium"],
			"Pa"  : [91, "Protactinium"],
			"U"   : [92, "Uranium"],
			"Np"  : [93, "Neptunium"],
			"Pu"  : [94, "Plutonium"],
			"Am"  : [95, "Americium"],
			"Cm"  : [96, "Curium"],
			"Bk"  : [97, "Berkelium"],
			"Cf"  : [98, "Californium"],
			"Es"  : [99, "Einsteinium"],
			"Fm"  : [100, "Fermium"],
			"Md"  : [101, "Mendelevium"],
			"No"  : [102, "Nobelium"],
			"Lr"  : [103, "Lawrencium"],
			"Rf"  : [104, "Rutherfordium"],
			"Db"  : [105, "Dubnium"],
			"Sg"  : [106, "Seaborgium"],
			"Bh"  : [107, "Bohrium"],
			"Hs"  : [108, "Hassium"],
			"Mt"  : [109, "Meitnerium"],
			"Ds"  : [110, "Darmstadtium"],
			"Rg"  : [111, "Roentgenium"],
			"Uub" : [112, "Ununbium"],
			"Uut" : [113, "Ununtrium"],
			"Uuq" : [114, "Ununquadium"],
			"Uup" : [115, "Ununpentium"],
			"Uuh" : [116, "Ununhexium"],
			"Uus" : [117, "Ununseptium"],
			"Uuo" : [118, "Ununoctium"]}

def isInt(string):
    try: 
        int(string)
        return True
    except ValueError:
        return False

def parseContracted(lines, i) :
	line = [s.strip() for s in lines[i].split()]
	primitives 	= int(line[0])
	orbitalType	= line[1]

	alpha 	= []
	c 		= []

	for j in range(i+1,i+primitives+1) :
		line = [s.strip() for s in lines[j].split()]
		alpha.append(float(line[0]))
		c.append(float(line[1]))

	returnStr = 'create_'+orbitalType.upper()+str(primitives)+'('
	for j in range(len(alpha)) :
		returnStr = returnStr+str(alpha[j])+','
	for j in range(len(c)-1) :
		returnStr = returnStr+str(c[j])+','
	return returnStr+str(c[-1])+');'


def parseName(line) :
	short = line[0]
	short = short.title()
	charge,name = names[short]
	basisName = line[1]
	return 'm_info = "'+name+' : '+basisName+'";'


if __name__ == '__main__' :
	if (not (len(sys.argv) == 2)) :
		print("Usage: python", __file__, "<basis file>")

	else :
		fileName = sys.argv[1]
		
		lines 				= []
		contracteds 		= []
		numberOfContracteds	= 0
		basisSize 			= -1
		firstBasisLine 		= -1
		nameLine 			= -1
		lastBasisLine 		= -1
		lineNumber			= -1

		with open(fileName, 'r') as inFile :
			lines = inFile.readlines()

		for i in range(len(lines)) :
			line = lines[i].strip().split()
			if line :
				if line[0] == "$basis" :
					firstBasisLine 	= i
					nameLine 		= firstBasisLine+2
					break

		nameString = parseName(lines[nameLine].split())

		for i in range(nameLine+1,len(lines)) :
			line = lines[i].split()
			if line :
				if isInt(line[0]) :
					if line[1].isalpha() :
						contracteds.append(parseContracted(lines,i))
						numberOfContracteds += 1

				elif line[0] == "$end" :
					break

		numberStr = 'setNumberOfOrbitals('+str(numberOfContracteds)+');'

		print(nameString)
		print(numberStr)
		for j in range(len(contracteds)) :
			print(contracteds[j])

