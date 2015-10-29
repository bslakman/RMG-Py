from rmgpy import settings
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.data.solvation import DatabaseError, SoluteData, SolvationDatabase, SolvationKinetics
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.reaction import *
from rmgpy.data.base import Database
#########################################################
database = Database()
species_dict = database.getSpecies('/Users/belinda/Code/RMG-Models/Liq_Fuel_Ox_2013/RMG_Dictionary.txt')

reaction_list = []

with open('/Users/belinda/Code/RMG-models/Liq_Fuel_Ox_2013/mechanism.txt', 'r') as mech_file:
    for line in mech_file:
        if line.strip().startswith('REACTIONS'): break
    for line in mech_file:
        if line.strip().startswith('!'): break
        if 'H_Abstraction' in line: reaction_list.append(line.strip())

kinetics_database = KineticsDatabase()
kinetics_database.load(os.path.join(settings['database.directory'], 'kinetics'), families=['H_Abstraction'], libraries=[])
family = kinetics_database.families['H_Abstraction']

for reaction in reaction_list:
    reactants, products = reaction.split('=')[0], reaction.split('=')[1]
    r1_ID, r2_ID = reactants.split('+')
    p1_ID, p2_ID = products.split('+')[0], products.split('+')[1]
    p2_ID = p2_ID.split(')')[0] + ')'
    r1 = species_dict[r1_ID]
    r2 = species_dict[r2_ID]
    p1 = species_dict[p1_ID]
    p2 = species_dict[p2_ID]
    reactant_molecules = [r1.molecule[0], r2.molecule[0]]
    product_molecules = [p1.molecule[0], p2.molecule[0]]
    testReaction = Reaction(reactants=[r1,r2], products=[p1,p2], reversible=True)

    checkRxn = kinetics_database.generateReactionsFromFamilies(reactant_molecules, product_molecules, only_families=[family])
    if checkRxn==[]:
    	reactIsoms = [mol.generateResonanceIsomers() for mol in reactant_molecules]
    	if len(reactIsoms)==2:
    		ri_1, ri_2 = reactIsoms
    		while checkRxn==[]:
    			for moleculeA in ri_1:
    				for moleculeB in ri_2:
    					reactant_molecules = [moleculeA, moleculeB]
    					checkRxn = kinetics_database.generateReactionsFromFamilies(reactant_molecules, product_molecules, only_families=[family])
#self.barrierDatabase = SolvationKinetics()
#        self.kineticsDatabase = KineticsDatabase()
#        self.kineticsDatabase.load(os.path.join(settings['database.directory'], 'kinetics'), families=['H_Abstraction'], libraries=[])
#        self.barrierDatabase.family = self.kineticsDatabase.families['H_Abstraction']
#        self.barrierDatabase.load(os.path.join(settings['database.directory'], 'kinetics', 'families', 'H_Abstraction'), None, None)
#
#
#        r1 = Species(molecule=[Molecule().fromAdjacencyList(
#"""
#1 *1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
#2 *2 H u0 p0 c0 {1,S}
#3    H u0 p0 c0 {1,S}
#4    H u0 p0 c0 {1,S}
#5    H u0 p0 c0 {1,S}
#"""
#        )]) # methane
#        r2 = Species(molecule=[Molecule().fromAdjacencyList(
#"""
#1 *3 O u1 p2 c0 {2,S}
#2    H u0 p0 c0 {1,S}
#"""
#        )]) # OH
#        p1 = Species(molecule=[Molecule().fromAdjacencyList(
#"""
#1 *3 C u1 p0 c0 {2,S} {3,S} {4,S}
#2    H u0 p0 c0 {1,S}
#3    H u0 p0 c0 {1,S}
#4    H u0 p0 c0 {1,S}
#"""
#        )]) # methyl
#        p2 = Species(molecule=[Molecule().fromAdjacencyList(
#"""
#1 *1 O u0 p2 c0 {2,S} {3,S}
#2 *2 H u0 p0 c0 {1,S}
#3    H u0 p0 c0 {1,S}
#"""
#        )]) # water
#        reaction = Reaction(reactants=[r1,r2], products=[p1,p2])
#        barrierCorrection = self.barrierDatabase.estimateBarrierCorrection(reaction)
