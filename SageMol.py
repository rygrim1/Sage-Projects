#!/usr/bin/python3
#Original code by Ryan Grimes written 3/15/20
from sage.plot.plot3d.shapes import LineSegment, Sphere
from sage.plot.colors import rainbow
from sage.plot.plot3d.implicit_plot3d import implicit_plot3d
from sage.plot.plot3d.shapes2 import text3d
from sage.repl.rich_output.pretty_print import show
from sage.calculus.var import var
from sage.plot.circle import circle
 
Radii = {'O': 0.73, 'N': 0.75, 'C': 0.77, 'He': 0.32, 'H': 0.37, 'S': 1.02, 'Cl': 0.99,

    'F': 0.71, 'Xe': 1.30, 'Si': 1.11, 'B': 0.82, 'P': 1.06, 'Br': 1.14}

Valency = {'O': 6, 'N': 5, 'C': 4, 'He': 0, 'H': 1, 'S': 6, 'Cl': 7, 

    'F': 7, 'Xe': 8, 'Si': 4, 'B': 3, 'P': 5, 'Br': 7}

Colors = {'O': 'red', 'N': 'darkblue', 'C': 'black', 'He': 'cyan', 'H': 'white', 'S': 'yellow',

	'Si': 'khaki', 'Xe': 'lightskyblue', 'F': 'green', 'Cl': 'limegreen', 'B': 'lightgoldenrodyellow', 'P': 'orange',

	'Br': 'darkred'}

Mol = Sphere(.00001, color = 'white').translate((0,0,0))

col = rainbow(10)

#-------------------------------------------------------------------------------------------------------------------------

def Planes(x,param,Key):

	yz = implicit_plot3d(lambda x,y,z: x, *param,color='white',opacity=0.8)

	xz = implicit_plot3d(lambda x,y,z: y, *param,color='cyan',opacity=0.8)

	xy = implicit_plot3d(lambda x,y,z: z, *param,color='grey',opacity=0.8)

	Plane = {'yz': yz, 'xz': xz, 'xy': xy}

	x += Plane[Key]

	return x

#-------------------------------------------------------------------------------------------------------------------------

def distance(A,B):

	D = ((B[0]-A[0])**2 + (B[1]-A[1])**2 + (B[2]-A[2])**2)**(1/2)

	return D

def Normalize(x):

	x = (x/0.77)*0.35

	return x

#-------------------------------------------------------------------------------------------------------------------------

def bent_geometry(XY2,Pi_Key,Pi_Bond):

	A = { 0 : (0,0,0) }

	S = { 2 : (0, 0.75, -0.75), 1 : (0, -0.75,-0.75) }

	C2R = { 1 : (0,0,1.25), 2 : (0,0,-2) }

	C2R_pair = { 1 : [C2R[1], C2R[2]] }

	C2R_label = { (0,0,1.5) : 'C2' }

	Central_Atom = { A[0] : 'S' }

	Central_Atom_Substitution = { A[0]: 'S' }

	Substituents = { S[1] : 'H', S[2] : 'H' }

	Additions = { S[1] : 'O', S[2] : 'O' }

	Sub_Sigma_Bonds = { 1 : [ A[0],S[1] ], 2 : [ A[0],S[2] ] }

	Sub_Pi_Bonds = { 1 : [ A[0],S[1] ], 2 : [ A[0],S[2] ] }

	Sym_Scale = { 'yz' : [(-1,1),(-1.25,1.25),(-1.5,0.75)], 'xz' : [(-1,1),(-1.25,1.25),(-1.5,0.75)] }

	Bent_Atoms = { }

	Axes = { }

	for sub in Additions:

		if sub in Substituents:

			del Substituents[sub]

	for sub in Central_Atom_Substitution:

		if sub in Central_Atom:

			del Central_Atom[sub]

	if Pi_Key == 'pi' and Pi_Bond == 'both':

		for sub in Sub_Pi_Bonds:

			if sub in Sub_Sigma_Bonds:

				del Sub_Sigma_Bonds[sub]

		for key, value in Sub_Pi_Bonds.items():

			XY2 += pi_planar(*value)

	if (Pi_Key == 'pi' and Pi_Bond == 1) or (Pi_Key == 'pi' and Pi_Bond ==2):

		del Sub_Sigma_Bonds[Pi_Bond]

		XY2 += pi_planar(*Sub_Pi_Bonds[Pi_Bond])

	Central_Atom.update(Central_Atom_Substitution)

	Substituents.update(Additions)

	if Substituents[S[1]] != Substituents[S[2]]:

		del (Sym_Scale['xz'], C2R_pair[1], C2R_label[(0,0,1.5)])

	Bent_Atoms = {**Substituents, **Central_Atom}
	
	for atom in Bent_Atoms:

		XY2 += Sphere(Normalize(Radii[Bent_Atoms[atom]]), color=Colors[Bent_Atoms[atom]]).translate(atom)

	for key, value in Sub_Sigma_Bonds.items():

		XY2 += LineSegment(*value,1, color = 'white')

	for key, value in C2R_pair.items():

		XY2 += LineSegment(*value,0.5, color = col.pop())

	for label in C2R_label:

		XY2 += text3d(C2R_label[label],label, color=col.pop())

	for plane in Sym_Scale:

		XY2 += Planes(XY2,Sym_Scale[plane],plane)

	return XY2

#-------------------------------------------------------------------------------------------------------------------------

def trig_planar(XY3):

	A = {0: (0,0,0)}

	S = {3: (0, 0.75, -0.70), 2: (0, -0.75,-0.70), 1: (0, 0, 1.026)}

	C3R = { 1 : (1.5, 0, 0), 2 : (-1.5, 0, 0) }

	C2R = { 1 : (0,0,1.75), 2 : (0,0,-1.5) }

	C3R_pairs = { 1: [C3R[1], C3R[2]] }

	C2R_pairs = { 2: [C2R[1], C2R[2]] }

	C2R_labels = { (0,0,1.825): 'C2' }

	C3R_labels = { (1.575,0,0): 'C3' }

	Substituents = { S[1]: 'H', S[2]: 'H', S[3]: 'H' }

	Additions = { S[1]: 'Cl', S[2]: 'H', S[3]: 'H'}

	Central_Atom = { A[0] : 'C' }

	Central_Atom_Substitution = { A[0]: 'B'}

	Sub_Bonds = { 'A0-S1' : [A[0],S[1]], 'A0-S2' : [A[0],S[2]], 'A0-S3' : [A[0],S[3]] }

	Sym_Planes = {1: 'yz', 2: 'xz'}

	Sym_Scale = {1: (-2,2), 2: (-2,2)}

	Trigonal_Planar_Atoms = { }

	Axes = { }

	Axes_labels = { }

	for sub in Additions:

		if sub in Substituents:

			del Substituents[sub]

	for sub in Central_Atom_Substitution:

		if sub in Central_Atom:

			del Central_Atom[sub]

	Central_Atom.update(Central_Atom_Substitution)

	Substituents.update(Additions)

	if Substituents[S[2]] != Substituents[S[3]]:

		del (Sym_Planes[2],C3R_pairs[1],C2R_pairs[2],C2R_labels[(0,0,1.825)],C3R_labels[(1.575,0,0)])

	Axes = {**C2R_pairs, **C3R_pairs}

	Axes_labels = {**C2R_labels, **C3R_labels}

	Trigonal_Planar_Atoms = {**Substituents, **Central_Atom}
	
	for atom in Trigonal_Planar_Atoms:

		XY3 += Sphere(Normalize(Radii[Trigonal_Planar_Atoms[atom]]), color=Colors[Trigonal_Planar_Atoms[atom]]).translate(atom)

	for key, value in Sub_Bonds.items():

		XY3 += LineSegment(*value,1, color = 'white')

	for key, value in Axes.items():

		XY3 += LineSegment(*value,0.5, color = col.pop())

	for label in Axes_labels:

		XY3 += text3d(Axes_labels[label],label, color=col.pop())

	Sym_Keys = list(Sym_Planes.values())

	XY3 += Planes(XY3,Sym_Scale[1],Sym_Keys)

	return XY3

#-------------------------------------------------------------------------------------------------------------------------

def tetra(XY4,q):

	#ADD PLANES OF SYMMETRY ALGORITHM

	A = { 0: (0,0,0) }

	S = { 3 : (2**(1/2), 1, 0), 2 : (-2**(1/2), 1, 0), 1: (0, -1, 2**(1/2)), 4: (0, -1, -2**(1/2)) }

	C2R = { 1 : (0,2.5,0), 2: (0,-2.5,0) }

	C3R = { 1 : (0, -2*.75, 2*2**(1/2)*.75), 2 : (0, 2*.75, -2*2**(1/2)*.75) }

	C2R_pair = { 1 : [C2R[1],C2R[2]] }

	C3R_pair = { 2 : [C3R[1],C3R[2]] }

	Axis_label = { C3R[1]: 'C3', (0,2.75,0): 'C2' }

	Substituents = { S[1]: 'H', S[2]: 'H', S[3]: 'H', S[4]: 'H' }

	Additions = { S[1]: 'H', S[2]: 'H', S[3]: 'H', S[4]: 'H' }

	Central_Atom = { A[0] : 'N' }

	Sub_Bonds = { 'A0-S1' : [A[0],S[1]], 'A0-S2' : [A[0],S[2]], 'A0-S3' : [A[0],S[3]], 'A0-S4': [A[0],S[4]] }

	Sym_Scale = { 'yz' : [(-1.5,1.5),(-2,2),(-2,2)], 'xy': [(-1.5,1.5),(-2,2),(-2,2)] }

	Axes = { }

	Tetra_Atoms = { }

	for sub in Additions:

		if sub in Substituents:

			del Substituents[sub]

	Substituents.update(Additions)

	if Substituents[S[3]] != Substituents[S[4]]:

		del (C3R_pair[2], Axis_label[C3R[1]])

	if Substituents[S[1]] != Substituents[S[4]]:

		del (C2R_pair[1], Axis_label[(0,2.75,0)], Sym_Scale['xy'])

	if q == 3:

		del (Substituents[S[1]], Sub_Bonds['A0-S1'], Axis_label[(0,2.5,0)], C2R_pair[1], Sym_Planes[2])

	if (Substituents[S[2]] != Substituents[S[3]]):

		del Sym_Scale['yz']

	Axes = {**C2R_pair, **C3R_pair}

	Tetra_Atoms = {**Substituents, **Central_Atom}
	
	for atom in Tetra_Atoms:

		XY4 += Sphere(Normalize(Radii[Tetra_Atoms[atom]]), color=Colors[Tetra_Atoms[atom]]).translate(atom)

	for key, value in Sub_Bonds.items():

		XY4 += LineSegment(*value,1, color = 'white')

	for key, value in Axes.items():	

		XY4 += LineSegment(*value,0.5, color = col.pop())

	for plane in Sym_Scale:

		XY4 += Planes(XY4,Sym_Scale[plane],plane)

	for label in Axis_label:

		XY4 += text3d(Axis_label[label],label, color=col.pop())

	return XY4

#-------------------------------------------------------------------------------------------------------------------------

def pi_ring_bond(Ax,Ay):

	delta = 0.05

	if Ax[1] == Ay[1]:
		
		Ax1 = (Ax[0], Ax[1]-delta, Ax[2])
		Ax2 = (Ax[0], Ax[1]+delta, Ax[2])
		Ay1 = (Ay[0], Ay[1]-delta, Ay[2])
		Ay2 = (Ay[0], Ay[1]+delta, Ay[2])

		PiBond = LineSegment(Ax1, Ay1, 1, color = 'white')
		PiBond += LineSegment(Ax2, Ay2, 1, color = 'white')

		return (PiBond)

	else:

		Ax1 = (Ax[0], Ax[1], Ax[2]-delta)
		Ax2 = (Ax[0], Ax[1], Ax[2]+delta)
		Ay1 = (Ay[0], Ay[1], Ay[2]-delta)
		Ay2 = (Ay[0], Ay[1], Ay[2]+delta)

		PiBond = LineSegment(Ax1, Ay1, 1, color = 'white')
		PiBond += LineSegment(Ax2, Ay2, 1, color = 'white')

		return (PiBond)

def pi_planar(Ax,Ay):

	if Ax[1] == Ay[1] and Ax[2] == Ay[2]:

		delta = 0.05

		Ax1 = (Ax[0], Ax[1], Ax[1]-delta)
		Ax2 = (Ax[0], Ax[1], Ax[1]+delta)
		Ay1 = (Ay[0], Ay[1], Ay[1]-delta)
		Ay2 = (Ay[0], Ay[1], Ay[1]+delta)

		PiBond = LineSegment(Ax1, Ay1, 1, color = 'white')
		PiBond += LineSegment(Ax2, Ay2, 1, color = 'white')

		return (PiBond)

	if Ay[1] == 0.75:

		delta = 0.05

		Ax1 = (Ax[0], Ax[1], Ax[1]-delta)
		Ax2 = (Ax[0], Ax[1], Ax[1]+delta)
		Ay1 = (Ay[0], Ay[1], -Ay[1]-delta)
		Ay2 = (Ay[0], Ay[1], -Ay[1]+delta)

		PiBond = LineSegment(Ax1, Ay1, 1, color = 'white')
		PiBond += LineSegment(Ax2, Ay2, 1, color = 'white')

		return (PiBond)

	if Ay[1] == -0.75:

		delta = 0.05

		Ax1 = (Ax[0], Ax[1]-delta, Ax[2])
		Ax2 = (Ax[0], Ax[1]+delta, Ax[2])
		Ay1 = (Ay[0], Ay[1]-delta, Ay[2])
		Ay2 = (Ay[0], Ay[1]+delta, Ay[2])

		PiBond = LineSegment(Ax1, Ay1, 1, color = 'white')
		PiBond += LineSegment(Ax2, Ay2, 1, color = 'white')

		return (PiBond)

#-------------------------------------------------------------------------------------------------------------------------

def planar_geometry(Mol,Key):

#-------------------------------------------------------------------------------------------------------------------------

	def Cyclobutadiene(Mol,A1,A2,A3,A4):

		Sub = { 1 : (1.50, 0, 0.75), 2 : (-1.50, 0, 0.75), 3 : (1.50, 0, -0.75-1.50), 4 : (-1.50, 0, -0.75-1.50) }

		C4R = { 1 : (0,2,-0.75), 2 : (0,-2,-0.75) }

		C4R_pair = { 1 : [C4R[1], C4R[2]] }

		C4R_label = { (0,2.25,-0.75) : 'C4'}

		Central_Bonds_Cyclo = { 'A1-A2' : [A1,A2], 'A1-A4' : [A1,A4], 'A2-A3' : [A2,A3], 'A3-A4' : [A3,A4] }

		Central_Atoms_Cyclo = { A1: 'C', A2: 'C', A3: 'C', A4: 'C' }

		Subs_Cyclo = { Sub[1] : 'H', Sub[2]: 'H', Sub[3]: 'H', Sub[4]: 'H'}

		Sub_Bonds_Cyclo = { 'A1-S1' : [A1,Sub[2]], 'A2-S2': [A2,Sub[1]], 'A3-S3': [A3,Sub[3]], 'A4-S4': [A4,Sub[4]] }

		Sym_Scale = { 'yz' : [(-1.5,1.5),(-2,2),(-2.5,1.0)], 'xz' : [(-2,2),(-2,2),(-2.5,1.0)] }

		for atom in Central_Atoms_Cyclo:

			Mol += Sphere(Normalize(Radii[Central_Atoms_Cyclo[atom]]), color=Colors[Central_Atoms_Cyclo[atom]]).translate(atom)

		for Substituent in Subs_Cyclo:

			Mol += Sphere(Normalize(Radii[Subs_Cyclo[Substituent]]), color=Colors[Subs_Cyclo[Substituent]]).translate(Substituent)

		for key, value in Sub_Bonds_Cyclo.items():	

			Mol += LineSegment(*value, 1, color='white')

		for key, value in Central_Bonds_Cyclo.items():

			Mol += LineSegment(*value, 2, color='white')

		for key, value in C4R_pair.items():

			Mol += LineSegment(*value,0.5, color = col.pop())

		for label in C4R_label:

			Mol += text3d(C4R_label[label],label, color=col.pop())

		for plane in Sym_Scale:

			Mol += Planes(Mol,Sym_Scale[plane],plane)

		#Mol += circle((1,0), 1, color='blue').plot3d()

		return Mol

#-------------------------------------------------------------------------------------------------------------------------

	def Ethylene(Mol,A1,A2):

		Sub = { 1 : (-1.50, 0, 0.75), 2 : (-1.50, 0, -0.75), 3: (1.50, 0, -0.75), 4: (1.50, 0, 0.75) }

		C2R = { 1 : (0,-2,0), 2 : (0,2,0), 3: (0,0,2), 4: (0,0,-2) }

		C2R_pair = { 1 : [C2R[1], C2R[2]], 2: [C2R[3],C2R[4]] }

		C2R_label = { (0,2.25,0) : 'C2', (0,0,2.25): 'C2'}

		Subs_Ethene = { Sub[1] : 'H', Sub[2]: 'Cl', Sub[3]: 'Cl', Sub[4]: 'H'}

		Central_Atoms_Ethene = { A1 : 'C', A2 : 'C' }

		Central_Bonds_Ethene = { 'A1-A2' : [A1,A2] } 

		Sub_Bonds_Ethene = { 'A1-S1' : [A1,Sub[1]], 'A1-S2' : [A1,Sub[2]], 'A2-S3' : [A2,Sub[3]], 'A2-S4': [A2,Sub[4]] }

		Sym_Scale = { 'yz' : [(-1.5,1.5),(-2,2),(-1,1)], 'xz' : [(-2,2),(-1,1),(-1,1)] }

		if Subs_Ethene[Sub[2]] == Subs_Ethene[Sub[3]]:

			del (C2R_pair[1], C2R_label[(0,2.25,0)])

		if (Subs_Ethene[Sub[2]] == Subs_Ethene[Sub[4]]) and (Subs_Ethene[Sub[1]] == Subs_Ethene[Sub[3]]):

			del (Sym_Scale['yz'], C2R_pair[2], C2R_label[(0,0,2.25)])

		Axes = {**C2R_pair}

		for atom in Central_Atoms_Ethene:

			Mol += Sphere(Normalize(Radii[Central_Atoms_Ethene[atom]]), color=Colors[Central_Atoms_Ethene[atom]]).translate(atom)

		for Substituent in Subs_Ethene:

			Mol += Sphere(Normalize(Radii[Subs_Ethene[Substituent]]), color=Colors[Subs_Ethene[Substituent]]).translate(Substituent)

		for key, value in Sub_Bonds_Ethene.items():	

			Mol += LineSegment(*value,1, color = 'white')

		for key, value in Central_Bonds_Ethene.items():	

			Mol += pi_planar(*value)

		for key, value in Axes.items():

			Mol += LineSegment(*value,0.5, color = col.pop())

		for label in C2R_label:

			Mol += text3d(C2R_label[label],label, color=col.pop())

		for plane in Sym_Scale:

			Mol += Planes(Mol,Sym_Scale[plane],plane)

		return Mol

	def Linear(Mol,A5,A0,A6):

		S = { 1 : (-2.25, 0, 0.75), 2 : (-2.25, 0, -0.75), 3: (2.25, -0.75, 0), 4: (2.25, 0.75, 0) }

		Linear_Subs = { S[1]: 'H', S[2]: 'H', S[3]: 'H', S[4]: 'H' }

		Linear_Atoms = { A5 : 'C', A0 : 'C', A6 : 'C' }

		Linear_Bonds = { 'A1-A0' : [A5,A0], 'A0-A2' : [A0,A6] }

		Sub_Bonds_Linear = { 'A1-S1': [A5,S[1]], 'A1-S2' : [A5,S[2]], 'A6-S3' : [A6,S[3]], 'A6-S4' : [A6,S[4]] }

		C2R = { 1 : (2.75,0,0), 2 : (-2.75,0,0), 3: (0,0,2), 4: (0,0,-2) }

		C2R_pair = { 1 : [C2R[1], C2R[2]], 2: [C2R[3],C2R[4]] }

		C2R_label = { (3,0,0) : 'C2', (0,0,2.25): 'C2'}

		Sym_Scale = { 'xz' : [(-2.50,2.50),(-1,1),(-1,1)] ,'xy' : [(-2.5,2.5),(-1,1),(-1,1)]}

		if (A5 != A0) or (A6 != A0):

			del (C2R_pair[2], C2R_label[(0,0,2.25)])

		Axes = {**C2R_pair}

		for Substituent in Linear_Subs:

			Mol += Sphere(Normalize(Radii[Linear_Subs[Substituent]]), color=Colors[Linear_Subs[Substituent]]).translate(Substituent)

		for atom in Linear_Atoms:

			Mol += Sphere(Normalize(Radii[Linear_Atoms[atom]]), color=Colors[Linear_Atoms[atom]]).translate(atom)

		for key, value in Linear_Bonds.items():	

			Mol += pi_planar(*value)

		for key, value in Sub_Bonds_Linear.items():	

			Mol += LineSegment(*value,1, color = 'white')

		for plane in Sym_Scale:

			Mol += Planes(Mol,Sym_Scale[plane],plane)

		for label in C2R_label:

			Mol += text3d(C2R_label[label],label, color=col.pop())

		for key, value in Axes.items():

			Mol += LineSegment(*value,0.5, color = col.pop())

		return Mol


	A = { 1: (-0.75, 0, 0), 2: (0.75, 0, 0), 4: (-0.75, 0, -1.50), 3: (0.75, 0, -1.50), 0: (0,0,0), 5: (-1.5,0,0), 6: (1.5,0,0) }

	if Key == 1:

		Mol += Ethylene(Mol,A[1],A[2])

	if Key == 2:

		Mol += Cyclobutadiene(Mol,A[1],A[2],A[3],A[4])

	if Key == 3:

		Mol += Linear(Mol,A[5],A[0],A[6])


	return Mol

#-------------------------------------------------------------------------------------------------------------------------

def Benzene():

	X = 0

	x = var('x')
	y = var('y')
	z = var('z')

	A = {1: (X, 0, 1.118), 2: (X, 1, 0.559), 3: (X, 1, -0.559), 4: (X, 0,-1.118), 5: (X, -1, -0.559), 6: (X, -1, 0.559) }

	S = {1: (1.75*X, 0, 1.75*1.118), 2: (1.75*X, 1.75*1, 1.75*0.559), 3: (1.75*X, 1.75*1, 1.75*-0.559), 4: (1.75*X, 0, 1.75*-1.118),

		5: (1.75*X, 1.75*-1, 1.75*-0.559), 6: (1.75*X, 1.75*-1, 1.75*0.559)}

	C2R = {1: (X, 0, 2), 2: (1.75*X, 1.25*1.75*1, 1.25*1.75*0.559), 3: (1.75*X, 1.25*1.75*1, 1.25*1.75*-0.559), 4: (X, 0, -2),

		5: (1.75*X, 1.25*1.75*-1, 1.25*1.75*-0.559), 6: (1.75*X, 1.25*1.75*-1, 1.25*1.75*0.559) }

	LC2 = { 1 : (X, 0, 2*1.1), 2 : (1.75*X, 1.1*1.25*1.75*1, 1.1*1.25*1.75*0.559), 3 : (1.75*X, 1.1*1.25*1.75*1, 1.1*1.25*1.75*-0.559), 

		4 : (X, 0, 1.1*-2), 5 : (1.75*X, 1.1*1.25*1.75*-1, 1.1*1.25*1.75*-0.559), 6 : (1.75*X, 1.1*1.25*1.75*-1, 1.1*1.25*1.75*0.559) }

	LC3 = {1: (X+2.25,0,0)}

	C3R = {1: (X-2,0,0), 2: (X+2,0,0)}

	O = (X,0,0)

	col = rainbow(10)

	Benzene = Sphere(.000001, color = 'white').translate(O)

	Sigma_Dict = { 'A1-A2' : [A[1],A[2]], 'A2-A3' : [A[2],A[3]], 'A3-A4' : [A[3],A[4]],

	'A4-A5' : [A[4],A[5]], 'A5-A6' : [A[5],A[6]], 'A1-A6' : [A[1],A[6]] }

	Pi_Dict = { 'A1-A2' : [A[1],A[2]], 'A5-A6' : [A[5],A[6]], 'A3-A4': [A[3],A[4]] }

	Atomic_Position = { A[1]: 'C', A[2]: 'C', A[3]: 'C', A[4]: 'C', A[5]: 'C', A[6]: 'C' }

	Sub_Position = { S[1]: 'H', S[2]: 'H', S[3]: 'H', S[4]: 'H', S[5]: 'H', S[6]: 'H' }

	Sub_Bonds = { 'A1-S1' : [A[1],S[1]], 'A2-S2' : [A[2],S[2]], 'A3-S3' : [A[3],S[3]], 'A4-S4' : [A[4],S[4]], 'A5-S5' : [A[5],S[5]],

		 'A6-S6' : [A[6],S[6]] }

	#Create common ring substituent dictionary ie OH, NO2 etc

	Subs = { S[2] : 'O', S[4]: 'H', S[6]: 'H'}

	#for i in range (0,len(C2R)):

	Axes_labels = { LC2[1]: 'C2', LC2[2]: 'C2', LC2[3]: 'C2', LC2[4]: 'C2', LC2[5]: 'C2', LC2[6]: 'C2' , LC3[1]: 'C3'}

	C2_Rotational_Pairs = { 'C2_1': [C2R[1],C2R[4]], 'C2_2': [C2R[2],C2R[5]], 'C2_3': [C2R[3], C2R[6]] }

	C3_Rotational_Pairs = { 'C6': [C3R[1], C3R[2]] }

	Symmetric_Planes = { 1 : 'yz', 2 : 'xz', 3 : 'xy'}

#-------------------------------------------------------------------------------------------------------------------------

	def Substituents(X,x):

		for substituent in Sub_Position:

			x += Sphere(Normalize(Radii[Sub_Position[substituent]]), color=Colors[Sub_Position[substituent]]).translate(substituent)

		for key, value in Sub_Bonds.items():

			x += LineSegment(*value,1, color = 'white')

		x += Sym(X,x)

		return x

	def Sym(X,x):

		x += Rotational_axes(X,x)
		x += Label_Baby_Jr(X,x)
		x += Planes(X,x)

		return x

	def Planes(X,x):

		col = ['white', 'cyan', 'aquamarine']

		x += Planar['yz']
		x += Planar['xz']

		#x += implicit_plot3d(lambda x,y,z: z, (-2,2), (-2,2), (-2,2),color=col.pop(),opacity=0.8)
		return x

	def Label_Baby_Jr(X,x):

		for label in Axes_labels:

			x += text3d(Axes_labels[label],label)

		return x

	def Rotational_axes(X,x):

		col = rainbow(8)

		for key, value in C2_Rotational_Pairs.items():	

			x += LineSegment(*value,0.5, color = 'blue')

		for key, value in C3_Rotational_Pairs.items():	

			x += LineSegment(*value,0.5, color = 'red')

		return x

	for atom in Atomic_Position:

		Benzene += Sphere(Normalize(Radii[Atomic_Position[atom]]), color=Colors[Atomic_Position[atom]]).translate(atom)

	for pi in Pi_Dict:

		if pi in Sigma_Dict:

			del Sigma_Dict[pi]

	for key, value in Sigma_Dict.items():

			Benzene += LineSegment(*value,1, color = 'white')

	for key, value in Pi_Dict.items():	

			Benzene += pi_ring_bond(*value)

	for sub in Subs:

		if sub in Sub_Position:

			del Sub_Position[sub]

	#Sigma_Dict.update(Pi_Dict)

	Sub_Position.update(Subs)

	show(Substituents(X,Benzene),frame=False)

#-------------------------------------------------------------------------------------------------------------------------

def Ethane(Mol,Key):

	x = 2**(1/2)

	A = { 1 : (0, -1, 0), 2 : (0, 1, 0) }

	S1 = { 1 : (-x, -2, -1), 2 : (x, -2, -1), 3: (0, -2, x ) }

	S2 = { 1 : (-x, 2, -1), 2 : (x, 2, -1), 3 : (0, 2, x ) }

	if Key == 'S':

		S2 = { 1 : (-x, 2, 1), 2 : (x, 2, 1), 3 : (0, 2, -x) }

	C2R = {1: [(0,0,1.75),(0,0,-1.75)], 2 : [(2,0,0),(-2,0,0)]}

	C3R = {'C3': [(0,-2.5,0),(0,2.5,0)]}

	Rot_label = { (0,0,2) : 'C2', (0,2.75,0) : 'C3', (2.25,0,0) : 'C2'}

	Substituents = { S1[1] : 'H', S1[2]: 'H', S1[3]: 'H', S2[1] : 'H', S2[2] : 'H', S2[3] : 'H' }

	Central_Atoms_Ethane = { A[1] : 'C', A[2] : 'C' }

	Central_Bonds_Ethane = { 'A1-A2' : [A[1],A[2]] } 

	Sub_Bonds_Ethane = { 'A1-S1' : [A[1],S1[1]], 'A1-S2' : [A[1],S1[2]], 'A1-S3' : [A[1],S1[3]], 'A2-S4': [A[2],S2[1]], 'A2-S5' : [A[2],S2[2]],
	
		'A2-S6': [A[2], S2[3]] }

	Sym_Scale = { 'yz' : [(-1.5,1.5),(-2,2),(-1,1)], 'xz': [(-1.5,1.5),(-2,2),(-1,1)] }

	if S2[1] != (-x, 2, -1):

		del (Sym_Scale['xz'], C2R[1], Rot_label[(0,0,2)])

	Atoms = {**Central_Atoms_Ethane, **Substituents} 

	Axes = {**C2R, **C3R}

	for atom in Atoms:

		Mol += Sphere(Normalize(Radii[Atoms[atom]]), color =Colors[Atoms[atom]]).translate(atom)

	Bonds = { **Central_Bonds_Ethane, **Sub_Bonds_Ethane }

	for key, value in Bonds.items():	

		Mol += LineSegment( *value, 1, color = 'white' )

	for plane in Sym_Scale:

		Mol += Planes(Mol,Sym_Scale[plane],plane)

	for key, value in Axes.items():	

		Mol += LineSegment(*value,0.5, color =col.pop())

	for label in Rot_label:

		Mol += text3d(Rot_label[label],label, color=col.pop())

	return Mol

#-------------------------------------------------------------------------------------------------------------------------

def Square_Pyramidal(Mol):

	A = { 0 : (0, 0 ,0) }

	S = { 1 : (0, 0, 1.5), 2 : (0, 1.5, 0), 3 : (0, -1.5, 0 ), 4 : (1.5,0,0), 5 : (-1.5,0,0)}

	C4R = {'C4': [(0,0,2.25),(0,0,-0.75)]}

	Axis_color = { 'C4' : 'blue', 'C2' : 'red' }

	Axis_label = { (0,0,2.50) : 'C4'}

	#C4R = {'C4': [(0,-1.5,0),(0,1.5,0)]}

	Substituents = { S[1] : 'F', S[2]: 'F', S[3]: 'F', S[4] : 'F', S[5] : 'F' }

	Central_Atom = { A[0] : 'Cl' }

	Bonds = { 'A0-S1' : [A[0],S[1]], 'A0-S2' : [A[0],S[2]], 'A0-S3' : [A[0],S[3]], 'A0-S4': [A[0],S[4]], 'A0-S5' : [A[0],S[5]] }

	Sym_Scale = { 'yz' : [ (-1.75,1.75),(-1.75,1.75),(-0.75,2) ], 'xz' : [ (-1.75,1.75),(-1.75,1.75),(-0.75,2) ] }

	Atoms = {**Central_Atom, **Substituents} 

	#Axes = {**C2R, **C4R }

	for atom in Atoms:

		Mol += Sphere(Normalize(Radii[Atoms[atom]]), color =Colors[Atoms[atom]]).translate(atom)

	for key, value in Bonds.items():	

		Mol += LineSegment( *value, 1, color = 'white' )

	for plane in Sym_Scale:

		Mol += Planes(Mol,Sym_Scale[plane],plane)

	for key, value in C4R.items():	

		Mol += LineSegment(*value,0.5, color = Axis_color[key])

	for label in Axis_label:

		Mol += text3d(Axis_label[label],label, color=col.pop())

	return Mol

#show(Square_Pyramidal(Mol))

#show(bent_geometry(Mol,'pi',1))

#

def trigonal_bipyramidal(Mol):

	A = {0: (0,0,0)}

	S = { 1: (0, 1.06, 0), 2: (0.75, -0.75, 0), 3: (-0.75,-0.75, 0), 4: (0,0, 1.06), 5: (0, 0, -1.06) }

	C3R = { 1 : (1.5, 0, 0), 2 : (-1.5, 0, 0) }

	C2R = { 1 : (0,0,1.75), 2 : (0,0,-1.5) }

	C3R_pairs = { 1: [C3R[1], C3R[2]] }

	C2R_pairs = { 2: [C2R[1], C2R[2]] }

	C2R_labels = { (0,0,1.825): 'C2' }

	C3R_labels = { (1.575,0,0): 'C3' }

	Substituents = { S[1] : 'Cl', S[2] : 'H', S[3] : 'H', S[4] : 'H', S[5] : 'H' }

	Additions = { S[1] : 'H', S[2] : 'H', S[3] : 'H', S[4] : 'H', S[5] : 'H' }

	Central_Atom = { A[0] : 'C' }

	Central_Atom_Substitution = { A[0]: 'B'}

	Bonds = { 'A0-S1' : [A[0],S[1]], 'A0-S2' : [A[0],S[2]], 'A0-S3' : [A[0],S[3]], 'A0-S4' : [A[0],S[4]], 'A0-S5': [A[0],S[5]] }

	Sym_Planes = { 1: 'yz' }

	Sym_Scale = { 'yz' : [(-1,1),(-1,1),(-1,1)] }

	Bipyramidal_Atoms = { **Substituents, **Central_Atom }

	Axes = { }

	Axes_labels = { }

	for atom in Bipyramidal_Atoms:

		Mol += Sphere(Normalize(Radii[Bipyramidal_Atoms[atom]]), color=Colors[Bipyramidal_Atoms[atom]]).translate(atom)

	for key, value in Bonds.items():	

		Mol += LineSegment( *value, 1, color = 'white' )

	for plane in Sym_Scale:

		Mol += Planes(Mol,Sym_Scale[plane],plane)

	for key, value in Axes.items():	

		Mol += LineSegment(*value,0.5, color =col.pop())

	return Mol

def Seesaw(Mol):

	A = { 0 : (0,0,0) }

	S = { 2 : (0.80, 0, -0.75), 1 : (-0.80, 0,-0.75), 3 : (0,-0.9,-0.25), 4 : (0,0.9,-0.25) }

	C2R = { 1 : (0,0,1.25), 2 : (0,0,-2) }

	C2R_pair = { 1 : [C2R[1], C2R[2]] }

	C2R_label = { (0,0,1.5) : 'C2' }

	Central_Atom = { A[0] : 'C' }

	Central_Atom_Substitution = { A[0]: 'S' }

	Substituents = { S[1] : 'H', S[2] : 'H', S[3] : 'He', S[4]: 'He'  }

	Additions = { S[1] : 'F', S[2] : 'F', S[3] : 'F', S[4]: 'F' }

	Bonds = { 1 : [ A[0],S[1] ], 2 : [ A[0],S[2] ], 3: [ A[0],S[3] ], 4: [A[0], S[4] ]  }

	Sym_Scale = { 'yz' : [(-1.25,1.25),(-1.25,1.25),(-1.25,0.75)], 'xz' : [(-1.25,1.25),(-1.25,1.25),(-1.25,0.75)] }

	Seesaw_Atoms = { }

	Axes = { }

	for sub in Additions:

		if sub in Substituents:

			del Substituents[sub]

	for sub in Central_Atom_Substitution:

		if sub in Central_Atom:

			del Central_Atom[sub]

	Central_Atom.update(Central_Atom_Substitution)

	Substituents.update(Additions)

	if Substituents[S[1]] != Substituents[S[2]]:

		del (Sym_Scale['xz'], C2R_pair[1], C2R_label[(0,0,1.5)])

	Atoms = {**Substituents, **Central_Atom}
	
	for atom in Atoms:

		Mol += Sphere(Normalize(Radii[Atoms[atom]]), color=Colors[Atoms[atom]]).translate(atom)

	for key, value in Bonds.items():

		Mol += LineSegment(*value,1, color = 'white')

	for key, value in C2R_pair.items():

		Mol += LineSegment(*value,0.5, color = col.pop())

	for label in C2R_label:

		Mol += text3d(C2R_label[label],label, color=col.pop())

	for plane in Sym_Scale:

		Mol += Planes(Mol,Sym_Scale[plane],plane)

	return Mol

def H2O2(Mol):

	x = 2**(1/2)

	A = { 1 : (0, -1, 0), 2 : (0, 1, 0) }

	S = { 1 : (-x, -2, -1), 2 : (x,2,1) }

	Axes = {'C2': [(0,0,1.75),(0,0,-1.75)]}

	C2R_label = { (0,0,2) : 'C2' }

	Substituents = { S[1] : 'H', S[2]: 'H' }

	Central_Atoms = { A[1] : 'O', A[2] : 'O' }

	Central_Bonds = { 'A1-A2' : [A[1],A[2]] } 

	Sub_Bonds = { 'A1-S1' : [A[1],S[1]], 'A2-S2' : [A[2],S[2]] }

	Atoms = { **Central_Atoms, **Substituents } 

	#for atom in Central_Atoms.values():
	#	if atom = A[1]

	if Central_Atoms[A[1]] == Central_Atoms[A[2]]:

		for key, value in Central_Bonds.items():

			Mol += LineSegment( *value, 1, color = Colors[Central_Atoms[A[1]]] )

	for atom in Atoms:

		Mol += Sphere(Normalize(Radii[Atoms[atom]]), color =Colors[Atoms[atom]]).translate(atom)

	#Bonds = { **Central_Bonds, **Sub_Bonds }

	for key, value in Sub_Bonds.items():	

		Mol += LineSegment( *value, 1, color = 'white' )

	for key, value in Axes.items():	

		Mol += LineSegment(*value,0.5, color =col.pop())

	for label in C2R_label:

		Mol += text3d(C2R_label[label],label, color=col.pop())

	return Mol


#show(planar_geometry(Mol,3),frame=False)

#show(H2O2(Mol), frame=True )

#show(tetra(Mol,1),frame=False)

#show(bent_geometry(Mol,'pi','both'),frame=False)

#show(Square_Pyramidal(Mol),frame=False)

#show(Ethane(Mol,'S'),frame=False)

#show(Seesaw(Mol),frame=False)

show(H2O2(Mol),frame=False)










