"""
Perform single point energy calculations in solvent on a set of gas phase
TS geometries in a database.
"""
import sqlite3
import openbabel as ob
from subprocess import Popen
import os
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
solvent = 'n-octane'

geometries = conn.execute("SELECT uniqueID, geometry FROM (SELECT * from TSData WHERE rxnFamily=\'" + family + "\') WHERE method='m062x'") # returns cursor
geometries = list(geometries.fetchall()) # make into list
entry = geometries[i-1]
#for n, entry in enumerate(geometries):
# Convert the geometry to xyz
obConversion = ob.OBConversion()
obConversion.SetInAndOutFormats("cml", "xyz")
mol = ob.OBMol()
obConversion.ReadString(mol, entry[1].encode('utf8'))
xyz_geom = obConversion.WriteString(mol)

if len(entry[0].encode('utf8').split('_')[0].split('+')) == 2:
	rmg_mol1 = Molecule().fromSMILES(entry[0].encode('utf8').split('_')[0].split('+')[0]) # for H-Abs, other bimolecular
	rmg_mol2 = Molecule().fromSMILES(entry[0].encode('utf8').split('_')[0].split('+')[1])
	mult = rmg_mol1.multiplicity + rmg_mol2.multiplicity - 1
else:
	rmg_mol = Molecule().fromSMILES(entry[0].encode('utf8').split('_')[0]) # for intra-H!
	mult = rmg_mol.multiplicity

xyz_geom = "0 {0}\n".format(mult) + '\n'.join(xyz_geom.split('\n')[2:])

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

# Check for successful output file; if not, do everything else
success = False
if os.path.exists(os.path.join(reaction_folder, "ts_" + solvent + ".log")):
	with open(os.path.join(reaction_folder, "ts_" + solvent + ".log")) as oF:
		for line in oF:
			line = line.strip()
			if "Normal termination" in line:
				success = True
				break
if not success:
	# Write input file
	options = "%mem=5GB\n%nprocshared=10"
	keywords = "# m062x/6-311+G(2df,2p) scrf(smd, solvent=" + solvent +") int=ultrafine freq nosymm"
	title = entry[0].encode('utf8')

	input_file_path = os.path.join(reaction_folder, "ts_" + solvent)
	input_file_ext = ".gjf"
	input_file = open(input_file_path + input_file_ext, 'w')
	input_file.write(options + "\n" + keywords + "\n\n" + title + "\n\n" + xyz_geom + "\n")
	input_file.close()

	# submits the input file to Gaussian
	#import ipdb; ipdb.set_trace()
	process = Popen(["/shared/apps/gaussian/G09_LINUX_LINDA/g09/g09", input_file_path + ".gjf", input_file_path + ".log"])
	process.communicate() # necessary to wait for executable termination!
