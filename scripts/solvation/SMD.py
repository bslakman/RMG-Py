#!/usr/bin/env python
# encoding: utf-8

# A script for gathering the SMD barrier correction values from Gaussian jobs.

import os.path
import argparse
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy import settings

parser = argparse.ArgumentParser(description="""
Given solvent, will get the delEa values from the output files.
""")
parser.add_argument("-s", "--solvent", help="Name of solvent, i.e. \'water\'")
parser.add_argument("-f", "--family", help="Name of family, i.e. \'intra_H_migration\'")
args= parser.parse_args()

rmgDatabase = RMGDatabase()
rmgDatabase.load(settings['database.directory'], kineticsFamilies=[args.family], reactionLibraries=[])

db_folder = '/Users/belinda/Code/RMG-database/input/kinetics/families'

data_dir = os.path.join('/Users/belinda/Gaussian/SMD/', args.family)
data_file = 'output_' + args.solvent + '.txt'
if not os.path.exists(os.path.join(db_folder, args.family, args.solvent, 'training')):
    os.makedirs(os.path.join(db_folder, args.family, args.solvent, 'training'))
output_file = os.path.join(db_folder, args.family, args.solvent, 'training', 'solvationReactions.py')
output_dict = os.path.join(db_folder, args.family, args.solvent, 'training', 'solvationDictionary.txt')

#print data_dir
#print data_file
#print output_file

out = open(output_file, 'w+')
out.write("""
#!/usr/bin/env python
# encoding: utf-8

name = \"""" + args.family + """\/training\"
shortDesc = u\"A list of solvation reactions used to train group additivity values\"
longDesc = u\"\"\"
\"\"\"
"""
)
out.close()

index = 1
unique = {}

for subdir, dirs, files in os.walk(data_dir):
    reactionString = os.path.basename(os.path.normpath(subdir))
    #if reactionString in ('H_Abs','solv2'): continue
    if not data_file in files: continue
    try:
        reactants, products = reactionString.split("_")
    except:
        print("String {0} is wrong length".format(reactionString))
        raise
    reactants = reactants.split("+")
    products = products.split("+")

    # Create the real reaction, so we can get the labeled atoms for the dictionary
    if args.family in ['H_Abstraction', 'Disproportionation', 'Cl-Abstraction']:
    	rSpecies1, rSpecies2 = [Species(molecule=[Molecule().fromSMILES(r)]) for r in reactants]
    	pSpecies1, pSpecies2 = [Species(molecule=[Molecule().fromSMILES(p)]) for p in products]
    	rSpecies1.generateResonanceIsomers()
    	rSpecies2.generateResonanceIsomers()
    	pSpecies1.generateResonanceIsomers()
    	pSpecies2.generateResonanceIsomers()
    	testReaction = Reaction(reactants=[rSpecies1, rSpecies2], products=[pSpecies1, pSpecies2], reversible=True)
    	reactionList = []
    	for moleculeA in rSpecies1.molecule:
    		for moleculeB in rSpecies2.molecule:
    			tempList = rmgDatabase.kinetics.generateReactionsFromFamilies([moleculeA, moleculeB], [], only_families=[args.family])
    			for rxn0 in tempList:
    				reactionList.append(rxn0)
    elif args.family in ['intra_H_migration']:
    	rSpecies = Species(molecule=[Molecule().fromSMILES(reactants[0])])
    	pSpecies = Species(molecule=[Molecule().fromSMILES(products[0])])
    	rSpecies.generateResonanceIsomers()
    	pSpecies.generateResonanceIsomers()
    	testReaction = Reaction(reactants=[rSpecies], products=[pSpecies], reversible=True)
    	reactionList = []
    	for moleculeA in rSpecies.molecule:
    		tempList = rmgDatabase.kinetics.generateReactionsFromFamilies([moleculeA], [], only_families=[args.family])
    		for rxn0 in tempList:
    			reactionList.append(rxn0)
    elif args.family in ['R_Addition_MultipleBond']:
        rSpecies1, rSpecies2 = [Species(molecule=Molecule().fromSMILES(r)) for r in reactants]
    	pSpecies = [Species(molecule=Molecule().fromSMILES(p)) for p in products]
    	rSpecies1.generateResonanceIsomers()
    	rSpecies2.generateResonanceIsomers()
    	pSpecies.generateResonanceIsomers()
    	testReaction = Reaction(reactants=[rSpecies1, rSpecies2], products=[pSpecies], reversible=False)
    	reactionList = []
    	for moleculeA in rSpecies1.molecule:
    		for moleculeB in rSpecies2.molecule:
    			tempList = rmgDatabase.kinetics.generateReactionsFromFamilies([moleculeA, moleculeB], [], only_families=[args.family])
    			for rxn0 in tempList:
    				reactionList.append(rxn0)
    gotOne=False
    for reaction in reactionList:
    	# Check if any of the RMG proposed reactions matches the reaction in the mechanism
    	if testReaction.isIsomorphic(reaction):
    		# Now add the labeled atoms to the Molecule, and check all labels were added
    		atLblsR = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])
    		atLblsP = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])

    		for reactant in reaction.reactants:
    			#reactant = reactant.molecule[0]
    			reactant.clearLabeledAtoms()
    			for atom in reactant.atoms:
    				for atomLabel in reaction.labeledAtoms:
    					if atom==atomLabel[1]:
    						atom.label = atomLabel[0]
    						atLblsR[atomLabel[0]] = True
    		for product in reaction.products:
    			#product = product.molecule[0]
    			product.clearLabeledAtoms()
    			for atom in product.atoms:
    				for atomLabel in reaction.labeledAtoms:
    					if atom==atomLabel[1]:
    						atom.label = atomLabel[0]
    						atLblsP[atomLabel[0]] = True
    		if all( atLblsR.values() ) and all( atLblsP.values() ):
    			gotOne=True
                my_reaction = reaction
                break

    if not gotOne: continue

    for i, r in enumerate(reactants):
        if r in unique:
            unique[r] += 1
            reactants[i] = "{0}-{1}".format(reactants[i], unique[r])
        else: unique.update({r: 1})
    for j, p in enumerate(products):
        if p in unique:
            unique[p] += 1
            products[j] = "{0}-{1}".format(products[j], unique[p])
        else: unique.update({p: 1})

    f = open(os.path.join(subdir,data_file), 'r')
    lines = f.readlines()
    f.close()
    for line in lines:
       if line.startswith("Difference"):
          text, barrierCorrection = line.split("= ")

    out = open(output_file, 'a')
    out.write('entry( \n')
    out.write('    index = ' + str(index) + ', \n')
    out.write('    reactants = [\'' + "\' , \'".join(reactants) + '\'], \n')
    out.write('    products = [\'' + "\' , \'".join(products) + '\'], \n')
    out.write('    solvent = \'' + args.solvent + '\', \n')
    out.write('    correction = BarrierCorrection(correction = (' + barrierCorrection + ', \'kJ/mol\')), \n')
    out.write('    shortDesc = u"""MO6-2X/MG3S calculations in g09 with SMD solvation model""", \n')
    out.write('    longDesc = \n u\"\"\" \n \"\"\" \n )')
    out.write('\n')
    out.close()

    index += 1

    # Create the dictionary
    if os.path.exists(output_dict):
        d = open(output_dict, 'r')
    else:
        d = open(output_dict, 'a+')
    entries = d.read().splitlines()
    d.close()
    for ind, reac in enumerate(reactants):
        if not reac in entries:
           d = open(output_dict, 'a')
           d.write(reac + '\n')
           d.write(my_reaction.reactants[ind].toAdjacencyList() + '\n')
           d.close()
    for indP, prod in enumerate(products):
        if not prod in entries:
            d = open(output_dict, 'a')
            d.write(prod + '\n')
            d.write(my_reaction.products[indP].toAdjacencyList() + '\n')
            d.close()
