import sys
import os
from rmgpy import settings
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.data.solvation import DatabaseError, SoluteData, SolvationDatabase, SolvationKinetics
from rmgpy.reaction import Reaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.data.base import Database
#########################################################
database = Database()
species_dict = database.getSpecies('/Users/belinda/Code/RMG-models/Liq_Fuel_Ox_2013/RMG_Dictionary.txt')

reaction_list = []

with open('/Users/belinda/Code/RMG-models/Liq_Fuel_Ox_2013/mechanism.txt', 'r') as mech_file:
    for line in mech_file:
        if line.strip().startswith('REACTIONS'): break
    for line in mech_file:
        if line.strip().startswith('!'): break
        if 'H_Abstraction' in line: reaction_list.append(line.strip())

rmg_database = RMGDatabase()
rmg_database.load(settings['database.directory'], kineticsFamilies=['H_Abstraction'], reactionLibraries=[])
kinetics_database = rmg_database.kinetics
family = kinetics_database.families['H_Abstraction']

barrier_database = SolvationKinetics()
barrier_database.family = family 
barrier_database.load(os.path.join(settings['database.directory'], 'kinetics', 'families', 'H_Abstraction'), None, None)

delEa_list = []

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

#reactant_molecules = [Molecule(SMILES="[CH3]"), Molecule(SMILES="O")]
#product_molecules = [Molecule(SMILES="C"), Molecule(SMILES="[OH]")] 
#
#reactant_species = [Species(molecule=mol.generateResonanceIsomers()) for mol in reactant_molecules]
#product_species = [Species(molecule=mol.generateResonanceIsomers()) for mol in product_molecules]
#
#testReaction = Reaction(reactants=reactant_species, products=product_species, reversible=True)

    reactions = kinetics_database.generateReactionsFromFamilies(reactant_molecules, product_molecules, only_families = [family.label])
    if len(reactions) == 0:
        reactIsoms = [mol.generateResonanceIsomers() for mol in reactant_molecules]
        if len(reactIsoms)==2:
            r1, r2 = reactIsoms
            checkRxn=[]
            #while checkRxn == []:
            for moleculeA in r1:
                for moleculeB in r2:
                    reactant_molecules = [moleculeA, moleculeB]
                    checkRxn = kinetics_database.generateReactionsFromFamilies(reactant_molecules, product_molecules, only_families=[family.label])
                    if checkRxn != []: break
            if checkRxn == []: 
                print testReaction.reactants, testReaction.products
                sys.exit()
            reaction = checkRxn[0]
    else:
        reaction = reactions[0]
    
    assert testReaction.isIsomorphic(reaction)
    
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
    
    
    barrierCorrection = barrier_database.estimateBarrierCorrection(reaction)
    delEa_list.append(barrierCorrection.correction.value_si) # Value in J/mol

new_mech_file = open('/Users/belinda/Code/RMG-models/Liq_Fuel_Ox_2013/mechanism_corrected.txt', 'a+')
num = 0 # Count H_Abs reactions
with open('/Users/belinda/Code/RMG-models/Liq_Fuel_Ox_2013/mechanism.txt', 'r') as mech_file:
    # write header and species list
    for line in mech_file:
        new_mech_file.write(line)
        if line.strip().startswith('REACTIONS'): break
    # write reactions, including the H_Abs one with corrected Ea
    for line in mech_file:
        if 'H_Abstraction' in line:
           # get the correction, apply to appropriate place in line
           corr = delEa_list[num]
           num += 1
           new_mech_file.write(line) 
        else: 
           new_mech_file.write(line)
           if line.strip().startswith('!'): break
    # write remainder of file
    for line in mech_file:
        new_mech_file.write(line)
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
