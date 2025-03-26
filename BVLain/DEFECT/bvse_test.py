from bvlain import Lain
import sys

file = sys.argv[1] # Give the structure as cif or vasp poscar file for bvse analysis
calc = Lain(verbose = False)
st = calc.read_file(file)

# Use this section for percolation radii analysis
#params = {'mobile_ion': 'Li1+',    # mobile specie
#		  'r_cut': 10.0,           # cutoff for interaction between the mobile species and framework
#		  'resolution': 0.2,	   # distance between the grid points
#}
#_ = calc.void_distribution(**params)
#radii = calc.percolation_radii()
#for key in radii.keys():
#    print(f'{key[-2:]} percolation radii is {round(radii[key], 4)} angstrom')

# Use this for percolation energy barrier analysis
params = {'mobile_ion': 'Li1+',    # mobile specie
		  'r_cut': 10.0,           # cutoff for interaction between the mobile species and framework
		  'resolution': 0.2,	   # distance between the grid points
		  'k': 100                 # maximum number of neighbors to be collected for each point
}
_ = calc.bvse_distribution(**params)
_ = calc.void_distribution(**params)

energies = calc.percolation_barriers(encut = 5.0)
for key in energies.keys():
    print(f'{key[-2:]} percolation barrier is {round(energies[key], 4)} eV')

calc.write_grd(file + '_bvse', task = 'bvse')  # saves .grd file
calc.write_cube(file + '_void', task = 'void') # save .cube file


# Use this for BV sum mismatch at every site in structure
#table = calc.mismatch(r_cut = 4.0)
#print(table.to_string())
