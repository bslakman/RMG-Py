"""
Perform single point energy calculations in solvent on a set of gas phase
reactant geometries (gets list from database and geometry from existing gas
phase calculations)
"""
import sqlite3
import openbabel as ob
from subprocess import Popen
import os
from shutil import copy
import sys
from rmgpy.molecule import Molecule

if len(sys.argv)>1:
	i = int(sys.argv[-1])
elif os.getenv('LSB_JOBINDEX'):
	i = int(os.getenv('LSB_JOBINDEX'))
else:
	raise Exception("Specify a TS number!")

conn = sqlite3.connect('/scratch/westgroup/ts_data.db') # Table name is TSData
print "Opened database successfully"

family = 'H_Abstraction'
solvent = 'water'

if family == 'H_Abstraction':
    geometries = conn.execute("SELECT uniqueID, r1, r2 FROM (SELECT * from TSData WHERE rxnFamily=\'" + family + "\') WHERE method='m062x'") # returns cursor
else:
	geometries = conn.execute("SELECT uniqueID, r1 FROM (SELECT * from TSData WHERE rxnFamily=\'" + family + "\') WHERE method='m062x'") # returns cursor
geometries = list(geometries.fetchall()) # make into list
entry = geometries[i-1]

# Check for folder
# Fine if reactants or products switched, too... for bimolecular
# It's not going to help if my SMILES are different from Pierre's though
reaction_folder = None
if len(entry[0].encode('utf8').split('_')[0].split('+')) == 2:
	r1 = entry[0].encode('utf8').split('_')[0].split('+')[0]
	r2 = entry[0].encode('utf8').split('_')[0].split('+')[1]
	p1 = entry[0].encode('utf8').split('_')[1].split('+')[0]
	p2 = entry[0].encode('utf8').split('_')[1].split('+')[1]
	possibleFolders = [os.path.join("/home/slakman.b/Gaussian/SMD/", family, "{0}+{1}_{2}+{3}".format(r1,r2,p1,p2)), os.path.join("/home/slakman.b/Gaussian/SMD/", family, "{0}+{1}_{2}+{3}".format(r1,r2,p2,p1)), os.path.join("/home/slakman.b/Gaussian/SMD/", family, "{0}+{1}_{2}+{3}".format(r2,r1,p1,p2)), os.path.join("/home/slakman.b/Gaussian/SMD/", family, "{0}+{1}_{2}+{3}".format(r2,r1,p2,p1))]
	for folder in possibleFolders:
	    if os.path.exists(folder):
			reaction_folder = folder
			break
else:
	folder = os.path.join("/home/slakman.b/Gaussian/SMD/", family, entry[0].encode('utf8'))
	if os.path.exists(folder): reaction_folder = folder

if reaction_folder is None:
	reaction_folder = os.path.join("/home/slakman.b/Gaussian/SMD/", family, entry[0].encode('utf8'))
	os.makedirs(reaction_folder)

err_message = None
# Find reactant gas-phase output file in scratch and copy to my reaction folder
if not os.path.exists(os.path.join(reaction_folder, entry[1] + ".log")):
    try:
    	gas_phase_output = os.path.join('/scratch/bhoorasingh.p/QMscratch/Species', entry[1], 'm062x', entry[1] + ".log")
    	copy(gas_phase_output, reaction_folder)
    except:
	    err_message = "Gas phase output file not found for {0}".format(entry[1])

if len(entry) > 2:
	if not os.path.exists(os.path.join(reaction_folder, entry[2] + ".log")):
		try:
			gas_phase_output2 = os.path.join('/scratch/bhoorasingh.p/QMscratch/Species', entry[2], 'm062x', entry[2] + ".log")
			copy(gas_phase_output2, reaction_folder)
		except:
			err_message = "Gas phase output file not found for {0}".format(entry[2])

gas_phase_output = os.path.join(reaction_folder, entry[1] + ".log")
gas_phase_output2 = None
if len(entry) > 2:
    gas_phase_output2 = os.path.join(reaction_folder, entry[2] + ".log")

# Find ts gas-phase output file in scratch and copy to my reaction folder
try:
	gas_phase_output_ts = os.path.join('/scratch/bhoorasingh.p/QMfiles/Reactions', family, entry[0], 'm062x', entry[0] + ".log")
	copy(gas_phase_output_ts, os.path.join(reaction_folder, "ts.log"))
except:
	err_message = "Gas phase output file not found for {0}".format(entry[0])

if err_message is None:
	if gas_phase_output2 is not None: gp_outputs = [[gas_phase_output, entry[1]], [gas_phase_output2, entry[2]]]
        else: gp_outputs = [[gas_phase_output, entry[1]]]
	# Extract geometry from gas-phase output
	for g in gp_outputs:
		gp = g[0]
		smiles = g[1]
		gotOutput = False
		# Check for liquid phase output in this or another folder.
		if os.path.exists(os.path.join(reaction_folder, smiles.encode('utf8') + "_" + solvent + ".log")):
			with open(os.path.join(reaction_folder, smiles.encode('utf8') + "_" + solvent + ".log")) as oF:
				for line in oF:
					line = line.strip()
					if "Normal termination" in line:
						gotOutput = True
					    break
			if gotOutput:
				continue

		for folder in os.listdir(os.path.join("/home/slakman.b/Gaussian/SMD/", family)):
			if folder == reaction_folder:
				continue
		    if os.path.exists(os.path.join("/home/slakman.b/Gaussian/SMD/", family, folder, smiles.encode('utf8') + "_" + solvent + ".log")):
		    	with open(os.path.join("/home/slakman.b/Gaussian/SMD/", family, folder, smiles.encode('utf8') + "_" + solvent + ".log")) as oF:
                	for line in oF:
                    	line = line.strip()
                        if "Normal termination" in line:
							copy(os.path.join("/home/slakman.b/Gaussian/SMD/", family, folder, smiles.encode('utf8') + "_" + solvent + ".log"), reaction_folder)
							gotOutput = True
							break
			if gotOutput=True:
				break

		if not gotOutput:
			xyz_geom = ""
			with open(gp, 'r') as gpf:
		            lines = gpf.read().split('\n')
		            input_geoms = [i for i,x in enumerate(lines) if 'Input orientation' in x]
			    geom_start = input_geoms[-1] + 5
			    geom_end = [lines[geom_start:].index(l) for l in lines[geom_start:] if "---" in l][0] + geom_start

			for atom_line in lines[geom_start:geom_end]:
			    atomic_number = int(atom_line.split()[1])
			    coordinates = atom_line.split()[-3:]
			    if atomic_number == 6: element = "C"
			    elif atomic_number == 1: element = "H"
			    elif atomic_number == 8: element = "O"
		            else: element = "x"
		            xyz_geom += "{0}\t{1} {2} {3}\n".format(element, coordinates[0], coordinates[1], coordinates[2])

			# multiplicity
			rmg_mol = Molecule().fromSMILES(smiles.encode('utf8'))
			mult = rmg_mol.multiplicity
			xyz_geom = "0 {0}\n".format(mult) + xyz_geom

			# Write solvation input file
			options = "%mem=5GB\n%nprocshared=10"
			keywords = "# m062x/6-311+G(2df,2p) scrf(smd, solvent=" + solvent +") int=ultrafine freq nosymm"
			title = smiles.encode('utf8')

			input_file_path = os.path.join(reaction_folder, smiles + "_" + solvent)
			input_file_ext = ".gjf"
			input_file = open(input_file_path + input_file_ext, 'w')
			input_file.write(options + "\n" + keywords + "\n\n" + title + "\n\n" + xyz_geom + "\n")
			input_file.close()

			# submits the input file to Gaussian
			process = Popen(["/shared/apps/gaussian/G09_LINUX_LINDA/g09/g09", input_file_path + ".gjf", input_file_path + ".log"])
			process.communicate() # necessary to wait for executable termination!

else:
    with open(os.path.join(reaction_folder, "err.txt"), 'w') as err_file:
        err_file.write(str(err_message))
