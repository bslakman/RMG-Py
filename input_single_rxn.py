import os
import sys

import rmgpy
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.data.kinetics import KineticsFamily, ReactionRecipe, KineticsDatabase
from rmgpy.data.rmg import RMGDatabase
from rmgpy.qm.main import QMCalculator
from rmgpy.qm.gaussian import GaussianTSB3LYP

if len(sys.argv)>1:
	i = int(sys.argv[-1])
elif os.getenv('LSB_JOBINDEX'):
	i = int(os.getenv('LSB_JOBINDEX'))
else:
	raise Exception("Specify a TS number!")

rxnFamiles = ['H_Abstraction']#,'Cl-Abstraction' ['intra_H_migration', 'R_Addition_MultipleBond', 'H_Abstraction', 'Disproportionation']

print 'Loading RMG Database ...'
rmgDatabase = RMGDatabase()
rmgDatabase.load(os.path.abspath(os.path.join(os.getenv('RMGpy'), '..', 'RMG-database', 'input')), kineticsFamilies=rxnFamiles)# solvation=False)
print 'Finished loading RMG Database ...'

print 'Loading solvation database'...
for family in rxnFamiles:
    barrier_database = SolvationKinetics()
	barrier_database.family = family
    barrier_database.load(os.path.join(os.getenv('RMGpy'), '..', 'RMG-database', 'input', 'kinetics', 'families', family), None, None)
print 'Finished loading solvation database for selected families'


rSpecies1 = Species(molecule=Molecule().fromSMILES("CCCCCCCC"))
rSpecies2 = Species(molecule=Molecule().fromSMILES("[O][O]"))
pSpecies1 = Species(molecule=Molecule().fromSMILES("[CH2]CCCCCCC"))
pSpecies2 = Species(molecule=Molecule().fromSMILES("O[O]"))
testReaction = Reaction(reactants=[rSpecies1, rSpecies2], products=[pSpecies1, pSpecies2], reversible=True)
reactionList = []
for moleculeA in rSpecies1.molecule:
	for moleculeB in rSpecies2.molecule:
		tempList = rmgDatabase.kinetics.generateReactionsFromFamilies([moleculeA, moleculeB], [], only_families=[rxnFamily])
		for rxn0 in tempList:
			reactionList.append(rxn0)

gotOne=False
for reaction in reactionList:
	# Check if any of the RMG proposed reactions matches the reaction in the mechanism
	if reaction.isIsomorphic(testReaction):
		# Now add the labeled atoms to the Molecule, and check all labels were added
		atLblsR = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])
		atLblsP = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])

		for reactant in reaction.reactants:
			reactant = reactant.molecule[0]
			reactant.clearLabeledAtoms()
			for atom in reactant.atoms:
				for atomLabel in reaction.labeledAtoms:
					if atom==atomLabel[1]:
						atom.label = atomLabel[0]
						atLblsR[atomLabel[0]] = True
		for product in reaction.products:
			product = product.molecule[0]
			product.clearLabeledAtoms()
			for atom in product.atoms:
				for atomLabel in reaction.labeledAtoms:
					if atom==atomLabel[1]:
						atom.label = atomLabel[0]
						atLblsP[atomLabel[0]] = True
		if all( atLblsR.values() ) and all( atLblsP.values() ):
			gotOne=True
			break

def calculate(reaction):
	rxnFamily = reaction.family.label
	tsDatabase = rmgDatabase.kinetics.families[rxnFamily].transitionStates
	# reaction = qmCalc.getKineticData(reaction, tsDatabase)
	# Need to get TS data, and if that works, get solvation data.
	solvationDatabase = rmgDatabase.kineticsfamilies[rxnFamily].solvationCorrections
	reaction = qmCalc.getSolvationData(reaction, tsDatabase, solvationDatabase)

	for files in os.listdir('./'):
		if files.startswith('core'):
			os.remove(files)

if not gotOne:
	print "No reactions found for reaction {4}: {0} + {1} = {2} + {3}".format(rSpecies1.molecule[0].toSMILES(), rSpecies2.molecule[0].toSMILES(), pSpecies1.molecule[0].toSMILES(), pSpecies2.molecule[0].toSMILES(), i)
else:
	qmCalc = QMCalculator(
									software='gaussian',
									method='m062x',
									fileStore='/scratch/slakman.b/QMfiles',
									scratchDirectory='/scratch/slakman.b/QMscratch',
									)
	calculate(reaction)
