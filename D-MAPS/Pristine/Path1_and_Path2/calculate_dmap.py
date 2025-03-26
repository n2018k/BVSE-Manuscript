# This is demo script to calculate D-Map for spodumene structure
# This example is for pristine spodumene structure with no defect along 16-20 path
# This is the lowest energy percolation energy path from BVLain
# Use Alpha_222.vasp as structure file, and bvparam file as parameter file
# The output from this is copied to plot.py for making the D-MAP plot in the manuscript

import sys, re, math
import numpy as np
from ase.io import read, write
from ase import Atom, Atoms
import matplotlib.pyplot as plt
import matplotlib.cm as cm

np.set_printoptions(threshold=10000000)

# read file with BV params
f = open(sys.argv[1], "r")
Lines = f.readlines()

# split lines out and store params
all_params = []
for line in Lines:
    temp = re.split('\s+', line)
    all_params.append(temp)

# Read the input structure
struct = read(sys.argv[2], format='vasp')

# get unique chemical symbols from struct
all_symbol = struct.get_chemical_symbols()
sym = [all_symbol[0]]

for s in all_symbol:
    if s not in sym:
        sym.append(s)

#ask valence for each unique atom for some reason here and store it
val = []

for s in sym:
    st = "Enter value for symbol: " + s + "\n"
    val.append(float(input(st)))


# define the 3d-grid for BV Sums
XStep = round(0.4/struct.get_cell_lengths_and_angles()[0],2)
YStep = round(1/struct.get_cell_lengths_and_angles()[1],2)
ZStep = round(1/struct.get_cell_lengths_and_angles()[2],2)

# now we need to make a map such that Li with correct required valence is used 
# to calculate bond valence sums

# for now its only for Li
# get BV params as average from the file

s_R0 = 0
s_b = 0
count = 0
allp = []
for i in all_params:
    if i[0] == 'Li' and i[2] == 'O':
        print(i)
        s_R0 += float(i[4]) # R0
        s_b += float(i[5]) # b
        count += 1 # count to get mean
        allp.append(i)

f_R0 = s_R0/count
f_b = s_b/count
print(f_R0, f_b)
# create the 3d grid

#Mapstep = [XStep, YStep, ZStep] # This wasnt used in this work
Mapstep = [0.02, 0.1, 0.02]
#Mapstep = [0.2, 0.2, 0.2]
print(Mapstep)
XGrid = [float(x/100) for x in range(0, 100, int(100*Mapstep[0]))]
YGrid = [float(x/100) for x in range(0, 1, int(100*Mapstep[1]))]
ZGrid = [float(x/100) for x in range(0, 101, int(100*Mapstep[2]))]

MapPts = []

for y in YGrid:
    for x in XGrid:
        for z in ZGrid:
            MapPts.append([x,y,z,0.0])

print(np.shape(MapPts))

# lets also make mappts translate as unitcell
# this will ease our life, we use that point as lithium atom and add it to box,
# calculate distances and bond valence sum

VMapPts = []

for m in MapPts:
    tt = 0
    dd = []
    new = struct.copy()
    Li_new=Atom('Li',[m[0]*struct.get_cell_lengths_and_angles()[0], 2.280092, m[2]*struct.get_cell_lengths_and_angles()[2]])
    new.append(Li_new)
    # calculate distance between this atom and all others 
    # and see if anything is too close
    for i in range(0, len(new)-1):
        d = new.get_distance(i,-1, mic=True)
        #print(new[i], i, d, new[-1])
        if new.get_chemical_symbols()[i] == 'O':
            dd.append(d)
        if d < 0.7:
            tt += 1
            break
    if tt == 0:
        # do the bond valence sum evaluation
        for d in dd:
            m[3] += np.exp((f_R0 - d) / f_b)
        VMapPts.append(m)
    

for m in MapPts:
    tt = 0
    dd = []
    new = struct.copy()
    Li_new=Atom('Li',[m[0]*struct.get_cell_lengths_and_angles()[0], 1.90199, m[2]*struct.get_cell_lengths_and_angles()[2]])
    new.append(Li_new)
    # calculate distance between this atom and all others 
    # and see if anything is too close
    for i in range(0, len(new)-1):
        d = new.get_distance(i,-1, mic=True)
        if new.get_chemical_symbols()[i] == 'O':
            dd.append(d)
        if d < 0.7:
            tt += 1
            break
    if tt == 0:
        # do the bond valence sum evaluation
        for d in dd:
            m[3] += np.exp((f_R0 - d) / f_b)
        VMapPts.append(m)

# I will get the bv sums out for each grid point here and use it separately to make plots

for n in range(1, 2):
	m4 = []
	m5 = []
	BVS = [] 
	temp = []
	for m in MapPts:
	   if m[1] == 0.0:
	       temp.append(m[3])
	       m4.append(m[0]*struct.get_cell_lengths_and_angles()[0])
	       m5.append(m[2]*struct.get_cell_lengths_and_angles()[2])
	       BVS.append(np.power(m[3]/2.43,-n))
        
	print("m4", m4)
	print("m5", m5)
	print("BVS", BVS)
	print(len(XGrid), len(ZGrid))

