#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""

"""
import rmgpy.quantity as quantity

import os
import os.path
import math
import logging
import numpy
from copy import copy, deepcopy

from base import *

import rmgpy.constants as constants
from rmgpy.molecule import Molecule, Atom, Bond, Group, atomTypes
from rmgpy.reaction import Reaction, ReactionError
from rmgpy.species import Species
from rmgpy.data.kinetics.common import UndeterminableKineticsError
################################################################################

def saveEntry(f, entry):
    """
    Write a Pythonic string representation of the given `entry` in the solvation
    database to the file object `f`.
    """
    f.write('entry(\n')
    f.write('    index = {0:d},\n'.format(entry.index))
    f.write('    label = "{0}",\n'.format(entry.label))
    
    if isinstance(entry.item, Molecule):
        if Molecule(SMILES=entry.item.toSMILES()).isIsomorphic(entry.item):
            # The SMILES representation accurately describes the molecule, so we can save it that way.
            f.write('    molecule = "{0}",\n'.format(entry.item.toSMILES()))
        else:
            f.write('    molecule = \n')
            f.write('"""\n')
            f.write(entry.item.toAdjacencyList(removeH=False))
            f.write('""",\n')
    elif isinstance(entry.item, Group):
        f.write('    group = \n')
        f.write('"""\n')
        f.write(entry.item.toAdjacencyList())
        f.write('""",\n')
    elif entry.item is not None:
        f.write('    group = "{0}",\n'.format(entry.item))
    
    if isinstance(entry.data, SoluteData):
        f.write('    solute = SoluteData(\n')
        f.write('        S = {0!r},\n'.format(entry.data.S))
        f.write('        B = {0!r},\n'.format(entry.data.B))
        f.write('        E = {0!r},\n'.format(entry.data.E))
        f.write('        L = {0!r},\n'.format(entry.data.L))
        f.write('        A = {0!r},\n'.format(entry.data.A))
        if entry.data.V is not None: f.write('        V = {0!r},\n'.format(entry.data.V))
        f.write('    ),\n')
    elif isinstance(entry.data, SolventData):
        f.write('    solvent = SolventData(\n')
        f.write('        s_g = {0!r},\n'.format(entry.data.s_g))
        f.write('        b_g = {0!r},\n'.format(entry.data.b_g))
        f.write('        e_g = {0!r},\n'.format(entry.data.e_g))
        f.write('        l_g = {0!r},\n'.format(entry.data.l_g))
        f.write('        a_g = {0!r},\n'.format(entry.data.a_g))
        f.write('        c_g = {0!r},\n'.format(entry.data.c_g))
        f.write('        s_h = {0!r},\n'.format(entry.data.s_h))
        f.write('        b_h = {0!r},\n'.format(entry.data.b_h))
        f.write('        e_h = {0!r},\n'.format(entry.data.e_h))
        f.write('        l_h = {0!r},\n'.format(entry.data.l_h))
        f.write('        a_h = {0!r},\n'.format(entry.data.a_h))
        f.write('        c_h = {0!r},\n'.format(entry.data.c_h))
        f.write('        A = {0!r},\n'.format(entry.data.A))
        f.write('        B = {0!r},\n'.format(entry.data.B))
        f.write('        C = {0!r},\n'.format(entry.data.C))
        f.write('        D = {0!r},\n'.format(entry.data.D))
        f.write('        E = {0!r},\n'.format(entry.data.E))
        f.write('        alpha = {0!r},\n'.format(entry.data.alpha))
        f.write('        beta = {0!r},\n'.format(entry.data.beta))
        f.write('        eps = {0!r},\n'.format(entry.data.eps))
        f.write('    ),\n')
    elif isinstance(entry.data, BarrierCorrection):
        f.write('    correction = BarrierCorrection(\n')
        if entry.data.correction is not None:
            f.write('        correction = ({0!r}, {1!r}),\n'.format(entry.data.correction.value_si, entry.data.correction.units))
        else:
            f.write('        correction = None,\n')
        f.write('    ),\n')
    elif entry.data is None:
        f.write('    solute = None,\n')
    else:
        raise DatabaseError("Not sure how to save {0!r}".format(entry.data))
    
    f.write('    shortDesc = u"""')
    try:
        f.write(entry.shortDesc.encode('utf-8'))
    except:
        f.write(entry.shortDesc.strip().encode('ascii', 'ignore')+ "\n")
    f.write('""",\n')
    f.write('    longDesc = \n')
    f.write('u"""\n')
    try:
        f.write(entry.longDesc.strip().encode('utf-8') + "\n")    
    except:
        f.write(entry.longDesc.strip().encode('ascii', 'ignore')+ "\n")
    f.write('""",\n')

    f.write(')\n\n')

def generateOldLibraryEntry(data):
    """
    Return a list of values used to save entries to the old-style RMG
    thermo database based on the thermodynamics object `data`.
    """
    raise NotImplementedError()
    
def processOldLibraryEntry(data):
    """
    Process a list of parameters `data` as read from an old-style RMG
    thermo database, returning the corresponding thermodynamics object.
    """
    raise NotImplementedError()

class SolvationReaction():
    """
    Stores information about kinetics of a solvation reaction
    """
    def __init__(self, reactants=None, products=None, solvent=None, gasKinetics=None, barrierCorrection=None):
        self.reactants = reactants #A list of species
	self.products = products #A list of species
        self.solvent = solvent #Solvent name 
        self.gasKinetics = gasKinetics #A kinetics object
        self.barrierCorrection = barrierCorrection #A quantity in energy units (kJ/mol, kcal/mol etc)

class BarrierCorrection():
    """
    Stores information about the correction to the energy between a reactant and its transition states.
    """
    def __init__(self, correction=None):
	self.correction = quantity.Energy(correction)
	self.comment = u''

class SolventData():
    """
    Stores Abraham/Mintz parameters for characterizing a solvent.
    Also associated with a category for solvation kinetics corrections. 
    """
    def __init__(self, s_h=None, b_h=None, e_h=None, l_h=None, a_h=None,
    c_h=None, s_g=None, b_g=None, e_g=None, l_g=None, a_g=None, c_g=None, A=None, B=None, 
    C=None, D=None, E=None, alpha=None, beta=None, eps=None, category=None):
        self.s_h = s_h
        self.b_h = b_h
        self.e_h = e_h
        self.l_h = l_h
        self.a_h = a_h
        self.c_h = c_h
        self.s_g = s_g
        self.b_g = b_g
        self.e_g = e_g
        self.l_g = l_g
        self.a_g = a_g
        self.c_g = c_g
        # These are parameters for calculating viscosity
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E
        # These are SOLUTE parameters used for intrinsic rate correction in H-abstraction rxns
        self.alpha = alpha
        self.beta = beta
        # This is the dielectric constant
        self.eps = eps
        # This is the 'category' of solvent it belogns to
        self.category = category
    
    def getHAbsCorrection(self):
        """
        If solvation is on, this will give the log10 of the ratio of the intrinsic rate
        constants log10(k_sol/k_gas) for H-abstraction rxns
        """
        return -8.3*self.alpha*self.beta
        
    def getSolventViscosity(self, T):
        """
        Returns the viscosity in Pa s, according to correlation in Perry's Handbook
        and coefficients in DIPPR
        """
        return math.exp(self.A + (self.B / T) + (self.C*math.log(T)) + (self.D * (T**self.E)))
                    
class SolvationCorrection():
    """
    Stores corrections for enthalpy, entropy, and Gibbs free energy when a species is solvated.
    Enthalpy and Gibbs free energy is in J/mol; entropy is in J/mol/K
    """
    def __init__(self, enthalpy=None, gibbs=None, entropy=None):
        self.enthalpy = enthalpy
        self.entropy = entropy
        self.gibbs = gibbs
            
class SoluteData():
    """
    Stores Abraham parameters to characterize a solute
    """
    def __init__(self, S=None, B=None, E=None, L=None, A=None, V=None, comment=""):
        self.S = S
        self.B = B
        self.E = E
        self.L = L
        self.A = A
        self.V = V
        self.comment = comment
    def __repr__(self):
        return "SoluteData(S={0},B={1},E={2},L={3},A={4},comment={5!r})".format(self.S, self.B, self.E, self.L, self.A, self.comment)
    
    def getStokesDiffusivity(self, T, solventViscosity):
        """
        Get diffusivity of solute using the Stokes-Einstein sphere relation. Radius is 
        found from the McGowan volume.
        """
        k_b = 1.3806488e-23 # m2*kg/s2/K
        radius = math.pow((75*self.V/3.14159/6.0221409e23),(1.0/3.0))/100 # in meters, V is in MgGowan volume in cm3/mol/100
        D = k_b*T/6/3.14159/solventViscosity/radius # m2/s
        return D
            
    def setMcGowanVolume(self, species):
        """
        Find and store the McGowan's Volume
        Returned volumes are in cm^3/mol/100 (see note below)
        See Table 2 in Abraham & McGowan, Chromatographia Vol. 23, No. 4, p. 243. April 1987
        doi: 10.1007/BF02311772
        
        "V is scaled to have similar values to the other
        descriptors by division by 100 and has units of (cm3mol−1/100)."
        the contibutions in this function are in cm3/mol, and the division by 100 is done at the very end.
        """
        molecule = species.molecule[0] # any will do, use the first.
        Vtot = 0

        for atom in molecule.atoms:
            thisV = 0.0
            if atom.isCarbon():
                thisV = 16.35
            elif (atom.element.number == 7): # nitrogen, do this way if we don't have an isElement method
                thisV = 14.39
            elif atom.isOxygen():
                thisV = 12.43
            elif atom.isHydrogen():
                thisV = 8.71
            elif (atom.element.number == 16):
                thisV = 22.91
            else:
                raise Exception()
            Vtot = Vtot + thisV

            for bond in molecule.getBonds(atom):
                # divide contribution in half since all bonds would be counted twice this way
                Vtot = Vtot - 6.56/2

        self.V= Vtot / 100; # division by 100 to get units correct.

################################################################################


################################################################################

class SolventLibrary(Database):
    """
    A class for working with a RMG solvent library.
    """
    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self,
                  index,
                  label,
                  solvent,
                  reference=None,
                  referenceType='',
                  shortDesc='',
                  longDesc='',
                  ):
        self.entries[label] = Entry(
            index = index,
            label = label,
            data = solvent,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
        )

    def load(self, path):
        """
        Load the solvent library from the given path
        """
        Database.load(self, path, local_context={'SolventData': SolventData}, global_context={})

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the solute database to the file object `f`.
        """
        return saveEntry(f, entry)
    
    def getSolventData(self, label):
        """
        Get a solvent's data from its name
        """
        return self.entries[label].data
        
        
class SoluteLibrary(Database):
    """
    A class for working with a RMG solute library
    """
    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self,
                  index,
                  label,
                  molecule,
                  solute,
                  reference=None,
                  referenceType='',
                  shortDesc='',
                  longDesc='',
                  ):
        try:
            mol = Molecule(SMILES=molecule)
        except:
            try:
                mol = Molecule().fromAdjacencyList(molecule)
            except:
                logging.error("Can't understand '{0}' in solute library '{1}'".format(molecule,self.name))
                raise

        self.entries[label] = Entry(
            index = index,
            label = label,
            item = mol,
            data = solute,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
        )
    
    def load(self, path):
        """
        Load the solute library from the given path
        """
        Database.load(self, path, local_context={'SoluteData': SoluteData}, global_context={})

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the solute database to the file object `f`.
        """
        return saveEntry(f, entry)

    def generateOldLibraryEntry(self, data):
        """
        Return a list of values used to save entries to the old-style RMG
        thermo database based on the thermodynamics object `data`.
        """
        return generateOldLibraryEntry(data)

    def processOldLibraryEntry(self, data):
        """
        Process a list of parameters `data` as read from an old-style RMG
        thermo database, returning the corresponding thermodynamics object.
        """
        return processOldLibraryEntry(data)

################################################################################

class SoluteGroups(Database):
    """
    A class for working with an RMG solute group additivity database.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def loadEntry(self,
                  index,
                  label,
                  group,
                  solute,
                  reference=None,
                  referenceType='',
                  shortDesc='',
                  longDesc='',
                  ):
        if group[0:3].upper() == 'OR{' or group[0:4].upper() == 'AND{' or group[0:7].upper() == 'NOT OR{' or group[0:8].upper() == 'NOT AND{':
            item = makeLogicNode(group)
        else:
            item = Group().fromAdjacencyList(group)
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = item,
            data = solute,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
        )
    
    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return saveEntry(f, entry)

    def generateOldLibraryEntry(self, data):
        """
        Return a list of values used to save entries to the old-style RMG
        thermo database based on the thermodynamics object `data`.
        """
        
        return generateOldLibraryEntry(data)

    def processOldLibraryEntry(self, data):
        """
        Process a list of parameters `data` as read from an old-style RMG
        thermo database, returning the corresponding thermodynamics object.
        """
        return processOldLibraryEntry(data)

################################################################################
class SolvationKinetics(Database):
    """
    A class for working with the RMG database of solvation kinetic corrections.
    """
    def __init__(self):
	self.groups = None
	self.family = None
	self.category = 'n-octane' 
    
    def load(self, path, local_context, global_context):
	if local_context is None: local_context = {}
        local_context['BarrierCorrection'] = BarrierCorrection

        fpath = os.path.join(path, self.category, 'training/solvationReactions.py')
        if os.path.exists(os.path.join(path, self.category)):
            depository = SolvationKineticsDepository(label='{0}/{1}/training'.format(path.split('/')[-1], self.category))
            depository.load(fpath, local_context, global_context )
            self.depository = depository
            
            fpath = os.path.join(path, self.category, 'solvationGroups.py')
            logging.debug("Loading solvation kinetics groups from {0}".format(fpath))
            groups = SolvationKineticsGroups(label='{0}/solvationGroups'.format(path.split('/')[-1]))
            groups.load(fpath, local_context, global_context)
            self.family.forwardTemplate.reactants = [groups.entries[entry.label] for entry in self.family.forwardTemplate.reactants]
            self.family.forwardTemplate.products = [groups.entries[entry.label] for entry in self.family.forwardTemplate.products]
            self.family.entries = groups.entries
            self.family.groups = groups
            groups.numReactants = len(self.family.forwardTemplate.reactants)
            self.groups = groups

    def estimateBarrierCorrection(self, reaction):
	return self.groups.estimateCorrectionUsingGroupAdditivity(reaction)

    def saveSolvationGroups(self, path, entryName='entry'):
        """
        Save the current database to the file at location `path` on disk. The
        optional `entryName` parameter specifies the identifier used for each
        data entry.
        """
        entries = self.groups.getEntriesToSave()

        # Write the header
        f = codecs.open(path, 'w', 'utf-8')
        f.write('#!/usr/bin/env python\n')
        f.write('# encoding: utf-8\n\n')
        f.write('name = "{0}"\n'.format(self.groups.name))
        f.write('shortDesc = u"{0}"\n'.format(self.groups.shortDesc))
        f.write('longDesc = u"""\n')
        f.write(self.groups.longDesc)
        f.write('\n"""\n\n')

        # Save the entries
        for entry in entries:
            saveEntry(f, entry)

        # Write the tree
        if len(self.groups.top) > 0:
            f.write('tree(\n')
            f.write('"""\n')
            f.write(self.generateOldTree(self.groups.top, 1))
            f.write('"""\n')
            f.write(')\n\n')

        f.close()

    def getForwardReactionForFamilyEntry(self, entry, family, groups, rxnFamily):
        """
        For a given `entry` for a reaction of the given reaction `family` (the
        string label of the family), return the reaction with the barrier correction
        for the "forward" direction as defined by the reaction family. For families 
        that are their own reverse, the direction the kinetics is given in will be 
        preserved. If the entry contains functional groups for the reactants, assume 
        that it is given in the forward direction and do nothing. Returns the reaction 
        in the direction consistent with the reaction family template, and the matching 
        template.
        """
        def matchSpeciesToMolecules(species, molecules):
            if len(species) == len(molecules) == 1:
                return species[0].isIsomorphic(molecules[0])
            elif len(species) == len(molecules) == 2:
                if species[0].isIsomorphic(molecules[0]) and species[1].isIsomorphic(molecules[1]):
                    return True
                elif species[0].isIsomorphic(molecules[1]) and species[1].isIsomorphic(molecules[0]):
                    return True
            return False
        
        reaction = None; template = None
        
        # Get the indicated reaction family
        if groups == None:
            raise ValueError('Invalid value "{0}" for family parameter.'.format(family))
            
        if all([(isinstance(reactant, Group) or isinstance(reactant, LogicNode)) for reactant in entry.item.reactants]):
            # The entry is a rate rule, containing functional groups only
            # By convention, these are always given in the forward direction and
            # have kinetics defined on a per-site basis
            reaction = Reaction(
                reactants = entry.item.reactants[:],
                products = [],
                kinetics = entry.data,
                degeneracy = 1,
            )
            template = [groups.entries[label] for label in entry.label.split(';')]
        
        elif (all([isinstance(reactant, (Molecule, Species)) for reactant in entry.item.reactants]) and
            all([isinstance(product, (Molecule, Species)) for product in entry.item.products])):
            # The entry is a real reaction, containing molecules
            # These could be defined for either the forward or reverse direction
            # and could have a reaction-path degeneracy
        
            reaction = Reaction(reactants=[], products=[])
            for molecule in entry.item.reactants:
                if isinstance(molecule, Molecule):
                    reactant = Species(molecule=[molecule])
                else:
                    reactant = molecule
                reactant.generateResonanceIsomers()
                reaction.reactants.append(reactant)
            for molecule in entry.item.products:
                if isinstance(molecule, Molecule):
                    product = Species(molecule=[molecule])
                else:
                    product = molecule
                product.generateResonanceIsomers()
                reaction.products.append(product)
            
            # Generate all possible reactions involving the reactant species
            generatedReactions = self.generateReactionsFromFamilies([reactant.molecule for reactant in reaction.reactants], [], only_families=[family], families=rxnFamily)
            # Remove from that set any reactions that don't produce the desired reactants and products
            forward = []; reverse = []
            for rxn in generatedReactions:
                if matchSpeciesToMolecules(reaction.reactants, rxn.reactants) and matchSpeciesToMolecules(reaction.products, rxn.products):
                    forward.append(rxn)
                if matchSpeciesToMolecules(reaction.reactants, rxn.products) and matchSpeciesToMolecules(reaction.products, rxn.reactants):
                    reverse.append(rxn)
            
            # We should now know whether the reaction is given in the forward or
            # reverse direction
            if len(forward) == 1 and len(reverse) == 0:
                # The reaction is in the forward direction, so use as-is
                reaction = forward[0]
                template = reaction.template
                # Don't forget to overwrite the estimated distances from the database with the distances for this entry
                reaction.correction = entry.data
            elif len(reverse) == 1 and len(forward) == 0:
                # The reaction is in the reverse direction
                # The reaction is in the forward direction, so use as-is
                reaction = reverse[0]
                template = reaction.template
                reaction.correction = entry.data
            elif len(reverse) > 0 and len(forward) > 0:
                print 'FAIL: Multiple reactions found for {0!r}.'.format(entry.label)
            elif len(reverse) == 0 and len(forward) == 0:
                print 'FAIL: No reactions found for "%s".' % (entry.label)
            else:
                print 'FAIL: Unable to estimate barrier correction for {0!r}.'.format(entry.label)
                
        assert reaction is not None
        assert template is not None
        return reaction, template
        
    def generateReactionsFromFamilies(self, reactants, products, only_families=None, families=None, **options):
        """
        Generate all reactions between the provided list of one or two
        `reactants`, which should be :class:`Molecule` objects. This method
        applies the reaction family.
        If `only_families` is a list of strings, only families with those labels
        are used.
        """
        # If there are two structures and they are the same, then make a copy
        # of the second one so we can independently manipulate both of them 
        # This is for the case where A + A --> products
        if len(reactants) == 2 and reactants[0] == reactants[1]:
            reactants[1] = reactants[1].copy(deep=True)
        
        reactionList = []
        reactionList.extend(families.generateReactions(reactants, **options))
        
        if products:
            reactionList = filterReactions(reactants, products, reactionList)
        
        return reactionList  
          
    def filterReactions(reactants, products, reactionList):
        """
        Remove any reactions from the given `reactionList` whose reactants do
        not involve all the given `reactants` or whose products do not involve 
        all the given `products`. This method checks both forward and reverse
        directions, and only filters out reactions that don't match either.
        """

        # Convert from molecules to species and generate resonance isomers.
        reactant_species = []
        for mol in reactants:
            s = Species(molecule=mol)
            s.generateResonanceIsomers()
            reactant_species.append(s)
        reactants = reactant_species
        product_species = []
        for mol in products:
            s = Species(molecule=mol)
            s.generateResonanceIsomers()
            product_species.append(s)
        products = product_species

        reactions = reactionList[:]

        for reaction in reactionList:
            # Forward direction
            reactants0 = [r for r in reaction.reactants]
            for reactant in reactants:
                for reactant0 in reactants0:
                    if reactant.isIsomorphic(reactant0):
                        reactants0.remove(reactant0)
                        break
            products0 = [p for p in reaction.products]
            for product in products:
                for product0 in products0:
                    if product.isIsomorphic(product0):
                        products0.remove(product0)
                        break
            forward = not (len(reactants0) != 0 or len(products0) != 0)
            # Reverse direction
            reactants0 = [r for r in reaction.products]
            for reactant in reactants:
                for reactant0 in reactants0:
                    if reactant.isIsomorphic(reactant0):
                        reactants0.remove(reactant0)
                        break
            products0 = [p for p in reaction.reactants]
            for product in products:
                for product0 in products0:
                    if product.isIsomorphic(product0):
                        products0.remove(product0)
                        break
            reverse = not (len(reactants0) != 0 or len(products0) != 0)
            if not forward and not reverse:
                reactions.remove(reaction)
        return reactions

class SolvationKineticsDepository(Database):
    
    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def __repr__(self):
        return '<SolvationKineticsDepository "{0}">'.format(self.label)
    
    def load(self, path, local_context=None, global_context=None):
        Database.load(self, path, local_context, global_context)
        
        # Load the species in the dictionary
        speciesDict = self.getSpecies(os.path.join(os.path.dirname(path),'solvationDictionary.txt'))
        # Make sure all of the reactions draw from only this set
        entries = self.entries.values()
        for entry in entries:
            # Create a new reaction per entry and convert the reactants and products to Species objects 
            # using the speciesDict
            rxn = entry.item 
            reactant_list = rxn.reactants
            product_list = rxn.products
            reactant_species = []
            product_species = []
            reversible = True
            assert reversible == rxn.reversible
            for reactant in reactant_list:
                if reactant not in speciesDict:
                    raise DatabaseError('Species {0} in kinetics depository {1} is missing from its dictionary.'.format(reactant, self.label))
                # For some reason we need molecule objects in the depository rather than species objects
                reactant_species.append(speciesDict[reactant])
            for product in product_list:
                product = product.strip()
                if product not in speciesDict:
                    raise DatabaseError('Species {0} in kinetics depository {1} is missing from its dictionary.'.format(product, self.label))
                # For some reason we need molecule objects in the depository rather than species objects
                product_species.append(speciesDict[product])
            rxn.reactants = reactant_species
            rxn.products = product_species
                
            if not rxn.isBalanced():
                raise DatabaseError('Reaction {0} in kinetics depository {1} was not balanced! Please reformulate.'.format(rxn, self.label))
                
    def loadEntry(self,
                  index,
                  reactants=None,
                  products=None,
                  solvent=None,
                  correction=None,
                  shortDesc='',
                  longDesc='',
                  ):
        reaction = Reaction(reactants=reactants, products=products, degeneracy=1, duplicate=False, reversible=True)
        entry = Entry(
            index = index,
            item = reaction,
            data = correction,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
        )
        self.entries[index] = entry
        return entry

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the database to the file object `f`.
        """
        return saveEntry(f, entry)

class SolvationKineticsGroups(Database):
    """
    A class for working with a database of group values for solvation kinetic corrections
    """

    def __init__(self, entries=None, top=None, label='', name='', shortDesc='', longDesc=''):
	Database.__init__(self, entries, top, label, name, shortDesc, longDesc)

    def __repr__(self):
	return '<SolvationKineticsGroups {0}:>'.format(self.label)

    def loadEntry(self, index, label, group, correction, reference=None, referenceType='', shortDesc='', longDesc=''):
        if group[0:3].upper() == 'OR{' or group[0:4].upper() == 'AND{' or group[0:7].upper() == 'NOT OR{' or group[0:8].upper() == 'NOT AND{':
            item = makeLogicNode(group)
        else:
            item = Group().fromAdjacencyList(group)
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = item,
            data = correction,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
        )

    def getReactionTemplate(self, reaction):
        """
        For a given `reaction` with properly-labeled :class:`Molecule` objects
        as the reactants, determine the most specific nodes in the tree that
        describe the reaction.
        """
        forwardTemplate = self.top[:]
        temporary = []
        symmetricTree = False
        for entry in forwardTemplate:
            if entry not in temporary:
                temporary.append(entry)
            else:
                # duplicate node found at top of tree
                # eg. R_recombination: ['Y_rad', 'Y_rad']
                assert len(forwardTemplate)==2 , 'Can currently only do symmetric trees with nothing else in them'
                symmetricTree = True
        forwardTemplate = temporary

        # Descend reactant trees as far as possible
        template = []
        for entry in forwardTemplate:
            # entry is a top-level node that should be matched
            group = entry.item

            # To sort out "union" groups, descend to the first child that's not a logical node
            # ...but this child may not match the structure.
            # eg. an R3 ring node will not match an R4 ring structure.
            # (but at least the first such child will contain fewest labels - we hope
            if isinstance(entry.item, LogicNode):
                group = entry.item.getPossibleStructures(self.entries)[0]

            atomList = group.getLabeledAtoms() # list of atom labels in highest non-union node

            for reactant in reaction.reactants:
                if isinstance(reactant, Species):
                    reactant = reactant.molecule[0]
                # Match labeled atoms
                # Check this reactant has each of the atom labels in this group
                if not all([reactant.containsLabeledAtom(label) for label in atomList]):
                    continue # don't try to match this structure - the atoms aren't there!
                # Match structures
                atoms = reactant.getLabeledAtoms()

                matched_node = self.descendTree(reactant, atoms, root=entry)
                if matched_node is not None:
                    template.append(matched_node)
                #else:
                #    logging.warning("Couldn't find match for {0} in {1}".format(entry,atomList))
                #    logging.warning(reactant.toAdjacencyList())

        # Get fresh templates (with duplicate nodes back in)
        forwardTemplate = self.top[:]
        if self.label.lower().startswith('r_recombination'):
            forwardTemplate.append(forwardTemplate[0])

        # Check that we were able to match the template.
        # template is a list of the actual matched nodes
        # forwardTemplate is a list of the top level nodes that should be matched
        if len(template) != len(forwardTemplate):
            logging.warning('Unable to find matching template for reaction {0} in reaction family {1}'.format(str(reaction), str(self)) )
            logging.warning(" Trying to match " + str(forwardTemplate))
            logging.warning(" Matched "+str(template))
            print str(self), template, forwardTemplate
            for n,reactant in enumerate(reaction.reactants):
                print "Reactant", n
                print reactant.toAdjacencyList() + '\n'
            for n,product in enumerate(reaction.products):
                print "Product", n
                print product.toAdjacencyList() + '\n'
            raise UndeterminableKineticsError(reaction)

        for reactant in reaction.reactants:
            if isinstance(reactant, Species):
                reactant = reactant.molecule[0]
            #reactant.clearLabeledAtoms()

        return template

    def estimateCorrectionUsingGroupAdditivity(self, reaction):
        """
        Determine the barrier correction due to solvation for a reaction 
        with the given `template` using group additivity.
        """
        template = self.getReactionTemplate(reaction)
        referenceCorrection = self.top[0].data # or something like that

        # Start with the generic barrier correction of the top-level nodes
        # Make a copy so we don't modify the original
        barrierCorrection = deepcopy(referenceCorrection)
        
        # Now add in more specific corrections if possible
        for node in template:
            entry = node
            comment_line = "Matched node "
            while not entry.data.correction and entry not in self.top:
                # Keep climbing tree until you find a (non-top) node with a correction 
                comment_line += "{0} >> ".format(entry.label)
                entry = entry.parent
            if entry.data.correction and entry not in self.top:
                barrierCorrection.correction.value_si += entry.data.correction.value_si
                comment_line += "{0} ({1})".format(entry.label, entry.longDesc.split('\n')[0])
            elif entry in self.top:
                comment_line += "{0} (Top node)".format(entry.label)
            barrierCorrection.comment += comment_line + '\n'
        
        return barrierCorrection
    
    def generateGroupAdditivityValues(self, trainingSet, user="Anonymous User"):
        """
        Generate the group additivity values using the given `trainingSet` of type (template, barrierCorrection). Returns 
	``True`` if the group values have changed significantly since the last time they
	 were fitted, or ``False`` otherwise.
        """
        # keep track of previous values so we can detect if they change
        old_entries = dict()
        for label,entry in self.entries.items():
            if entry.data is not None:
                old_entries[label] = entry.data

        # Determine a complete list of the entries in the database, sorted as in the tree
        groupEntries = self.top[:]

        for entry in self.top:
            groupEntries.extend(self.descendants(entry)) # Entries in the solvation kinetics groups tree

        # Determine a unique list of the groups we will be able to fit parameters for
        groupList = []
        for template, corrections in trainingSet:
            for group in template:
                if group not in self.top:
                    groupList.append(group)
                    groupList.extend(self.ancestors(group)[:-1])
        groupList = list(set(groupList))
        groupList.sort(key=lambda x: x.index)

        if True: # should remove this IF block, as we only have one method.
            # Initialize dictionaries of fitted group values and uncertainties
            groupValues = {}; groupUncertainties = {}; groupCounts = {}; groupComments = {}
            for entry in groupEntries:
                groupValues[entry] = None 
                groupUncertainties[entry] = []
                groupCounts[entry] = []
                groupComments[entry.label] = set()

            # Generate least-squares matrix and vector
            A = []; b = []

            correction_data = []
	    for template, correctionData in trainingSet:
                d = correctionData.correction.value_si
                correction_data.append(d)

                # Create every combination of each group and its ancestors with each other
                combinations = []
                for group in template:
                    groups = [group]; groups.extend(self.ancestors(group)) # Groups from the group.py tree
                    combinations.append(groups)
                combinations = getAllCombinations(combinations)
                # Add a row to the matrix for each combination
                for groups in combinations:
                    Arow = [1 if group in groups else 0 for group in groupList]
                    Arow.append(1)
                    brow = d
                    A.append(Arow); b.append(brow)

                    for group in groups:
                        groupComments[group.label].add("{0!s}".format(template))

            if len(A) == 0:
                logging.warning('Unable to fit kinetics groups for family "{0}"; no valid data found.'.format(self.label))
                return
            A = numpy.array(A)
            b = numpy.array(b)
            correction_data = numpy.array(correction_data)

            x, residues, rank, s = numpy.linalg.lstsq(A, b)
            
            # Determine error in each group
            stdev = numpy.zeros(len(groupList)+1, numpy.float64)
            count = numpy.zeros(len(groupList)+1, numpy.int)

            for index in range(len(trainingSet)):
                template, correction = trainingSet[index]
                d = numpy.float64(correction_data[index])
                dm = x[-1] + sum([x[groupList.index(group)] for group in template if group in groupList])
                variance = (dm - d)**2
                for group in template:
                    groups = [group]
                    groups.extend(self.ancestors(group))
                    for g in groups:
                        if g.label not in [top.label for top in self.top]:
                            ind = groupList.index(g)
                            stdev[ind] += variance
                            count[ind] += 1
                stdev[-1] += variance
                count[-1] += 1

            import scipy.stats
            ci = numpy.zeros(len(count))
            for i in range(len(count)):
                if count[i] > 1:
                    stdev[i] = numpy.sqrt(stdev[i] / (count[i] - 1))
                    ci[i] = scipy.stats.t.ppf(0.975, count[i] - 1) * stdev[i]
                else:
                    stdev[i] = None
                    ci[i] = None
            # Update dictionaries of fitted group values and uncertainties
            for entry in groupEntries:
                if entry == self.top[0]:
                    groupValues[entry]=x[-1]
                    groupUncertainties[entry].append(ci[-1])
                    groupCounts[entry].append(count[-1])
                elif entry.label in [group.label for group in groupList]:
                    index = [group.label for group in groupList].index(entry.label)
                    groupValues[entry]=x[index]
                    groupUncertainties[entry].append(ci[index])
                    groupCounts[entry].append(count[index])
                else:
                    groupValues[entry] = None
                    groupUncertainties[entry] = None
                    groupCounts[entry] = None

            # Store the fitted group values and uncertainties on the associated entries
            for entry in groupEntries:
                if groupValues[entry] is not None:
                    if not any(numpy.isnan(numpy.array(groupUncertainties[entry]))):
                        # should be entry.data.* (e.g. entry.data.uncertainties)
                        uncertainties = numpy.array(groupUncertainties[entry])
                        uncertaintyType = '+|-'
                    else:
                        uncertainties = {}
                    # should be entry.*
                    shortDesc = "Fitted to {0} solvation corrections.\n".format(groupCounts[entry][0])
		    longDesc = "\n".join(groupComments[entry.label])
                    # we don't include uncertainties for the moment (not in BarrierCorrection object)
                    entry.data = BarrierCorrection(correction= (groupValues[entry], 'J/mol'))
                    entry.shortDesc = shortDesc
                    entry.longDesc = longDesc
                else:
                    entry.data = BarrierCorrection()

        # Add a note to the history of each changed item indicating that we've generated new group values
        import time
        event = [time.asctime(), user, 'action', 'Generated new group additivity values for this entry.']
        changed = False
        for label, entry in self.entries.items():
            if entry.data is not None:
                continue # because this is broken:
                if old_entries.has_key(label):
                    old_entry = old_entries[label][label][0]
                    for key, distance in entry.data.iteritems():
                        diff = 0
                        for k in range(3):
                            diff += abs(distance[0][k]/old_entry[k] - 1)
                        if diff > 0.01:
                            changed = True
                            entry.history.append(event)
            else:
                changed = True
                entry.history.append(event)
        return True # because the thing above is broken
        return changed

################################################################################
class SolvationDatabase(object):
    """
    A class for working with the RMG solvation thermodynamics database.
    """

    def __init__(self):
        self.libraries = {}
        self.libraries['solvent'] = SolventLibrary()
        self.libraries['solute'] = SoluteLibrary()
        self.groups = {}
        self.local_context = {
            'SoluteData': SoluteData,
            'SolventData': SolventData
        }
        self.global_context = {}

    def __reduce__(self):
        """
        A helper function used when pickling a SolvationDatabase object.
        """
        d = {
            'libraries': self.libraries,
            'groups': self.groups,
            }
        return (SolvationDatabase, (), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling a SolvationDatabase object.
        """
        self.libraries = d['libraries']
        self.groups = d['groups']

    def load(self, path, libraries=None, depository=True):
        """
        Load the solvation database from the given `path` on disk, where `path`
        points to the top-level folder of the solvation database.
        
        Load the solvent and solute libraries, then the solute groups.
        """
        
        self.libraries['solvent'].load(os.path.join(path,'libraries','solvent.py'))
        self.libraries['solute'].load(os.path.join(path,'libraries','solute.py'))
         
        self.loadGroups(os.path.join(path, 'groups'))
        
    def getSolventData(self, solvent_name):
        try:
            solventData = self.libraries['solvent'].getSolventData(solvent_name)
        except:
            raise DatabaseError('Solvent {0!r} not found in database'.format(solvent_name))
        return solventData
        
        
    def loadGroups(self, path):
        """
        Load the solute database from the given `path` on disk, where `path`
        points to the top-level folder of the solute database.
        
        Three sets of groups for additivity, atom-centered ('abraham'), non atom-centered 
        ('nonacentered'), and radical corrections ('radical')
        """
        logging.info('Loading Platts additivity group database from {0}...'.format(path))
        self.groups = {}
        self.groups['abraham']   =   SoluteGroups(label='abraham').load(os.path.join(path, 'abraham.py'  ), self.local_context, self.global_context)
        self.groups['nonacentered']  =  SoluteGroups(label='nonacentered').load(os.path.join(path, 'nonacentered.py' ), self.local_context, self.global_context)
        self.groups['radical']  =  SoluteGroups(label='radical').load(os.path.join(path, 'radical.py' ), self.local_context, self.global_context)
   
    def save(self, path):
        """
        Save the solvation database to the given `path` on disk, where `path`
        points to the top-level folder of the solvation database.
        """
        path = os.path.abspath(path)
        if not os.path.exists(path): os.mkdir(path)
        self.saveLibraries(os.path.join(path, 'libraries'))
        self.saveGroups(os.path.join(path, 'groups'))

    def saveLibraries(self, path):
        """
        Save the solute libraries to the given `path` on disk, where `path`
        points to the top-level folder of the solute libraries.
        """
        if not os.path.exists(path): os.mkdir(path)
        for library in self.libraries.keys():
            self.libraries[library].save(os.path.join(path, library+'.py'))
        
    def saveGroups(self, path):
        """
        Save the solute groups to the given `path` on disk, where `path`
        points to the top-level folder of the solute groups.
        """
        if not os.path.exists(path): os.mkdir(path)
        for group in self.groups.keys():
            self.groups[group].save(os.path.join(path, group+'.py'))

    def loadOld(self, path):
        """
        Load the old RMG solute database from the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """
        
        for (root, dirs, files) in os.walk(os.path.join(path, 'thermo_libraries')):
            if os.path.exists(os.path.join(root, 'Dictionary.txt')) and os.path.exists(os.path.join(root, 'Library.txt')):
                library = SoluteLibrary(label=os.path.basename(root), name=os.path.basename(root))
                library.loadOld(
                    dictstr = os.path.join(root, 'Dictionary.txt'),
                    treestr = '',
                    libstr = os.path.join(root, 'Library.txt'),
                    numParameters = 5,
                    numLabels = 1,
                    pattern = False,
                )
                library.label = os.path.basename(root)
                self.libraries[library.label] = library

        self.groups = {}
        self.groups['abraham'] = SoluteGroups(label='abraham', name='Platts Group Additivity Values for Abraham Solute Descriptors').loadOld(
            dictstr = os.path.join(path, 'thermo_groups', 'Abraham_Dictionary.txt'),
            treestr = os.path.join(path, 'thermo_groups', 'Abraham_Tree.txt'),
            libstr = os.path.join(path, 'thermo_groups', 'Abraham_Library.txt'),
            numParameters = 5,
            numLabels = 1,
            pattern = True,
        )

    def saveOld(self, path):
        """
        Save the old RMG Abraham database to the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """
        # Depository not used in old database, so it is not saved

        librariesPath = os.path.join(path, 'thermo_libraries')
        if not os.path.exists(librariesPath): os.mkdir(librariesPath)
        for library in self.libraries.values():
            libraryPath = os.path.join(librariesPath, library.label)
            if not os.path.exists(libraryPath): os.mkdir(libraryPath)
            library.saveOld(
                dictstr = os.path.join(libraryPath, 'Dictionary.txt'),
                treestr = '',
                libstr = os.path.join(libraryPath, 'Library.txt'),
            )

        groupsPath = os.path.join(path, 'thermo_groups')
        if not os.path.exists(groupsPath): os.mkdir(groupsPath)
        self.groups['abraham'].saveOld(
            dictstr = os.path.join(groupsPath, 'Abraham_Dictionary.txt'),
            treestr = os.path.join(groupsPath, 'Abraham_Tree.txt'),
            libstr = os.path.join(groupsPath, 'Abraham_Library.txt'),
        )

    def getSoluteData(self, species):
        """
        Return the solute descriptors for a given :class:`Species`
        object `species`. This function first searches the loaded libraries
        in order, returning the first match found, before falling back to
        estimation via Platts group additivity.
        """
        soluteData = None
        
        # Check the library first
        soluteData = self.getSoluteDataFromLibrary(species, self.libraries['solute'])
        if soluteData is not None:
            assert len(soluteData)==3, "soluteData should be a tuple (soluteData, library, entry)"
            soluteData[0].comment += "Data from solute library"
            soluteData = soluteData[0]
        else:
            # Solute not found in any loaded libraries, so estimate
            soluteData = self.getSoluteDataFromGroups(species)
            # No Platts group additivity for V, so set using atom sizes
            soluteData.setMcGowanVolume(species)
        # Return the resulting solute parameters S, B, E, L, A
        return soluteData

    def getAllSoluteData(self, species):
        """
        Return all possible sets of Abraham solute descriptors for a given
        :class:`Species` object `species`. The hits from the library come
        first, then the group additivity  estimate. This method is useful
        for a generic search job. Right now, there should either be 1 or 
        2 sets of descriptors, depending on whether or not we have a 
        library entry.
        """
        soluteDataList = []
        
        # Data from solute library
        data = self.getSoluteDataFromLibrary(species, self.libraries['solute'])
        if data is not None: 
            assert len(data) == 3, "soluteData should be a tuple (soluteData, library, entry)"
            data[0].comment += "Data from solute library"
            soluteDataList.append(data)
        # Estimate from group additivity
        # Make it a tuple
        data = (self.getSoluteDataFromGroups(species), None, None)
        soluteDataList.append(data)
        return soluteDataList

    def getSoluteDataFromLibrary(self, species, library):
        """
        Return the set of Abraham solute descriptors corresponding to a given
        :class:`Species` object `species` from the specified solute
        `library`. If `library` is a string, the list of libraries is searched
        for a library with that name. If no match is found in that library,
        ``None`` is returned. If no corresponding library is found, a
        :class:`DatabaseError` is raised.
        """
        for label, entry in library.entries.iteritems():
            for molecule in species.molecule:
                if molecule.isIsomorphic(entry.item) and entry.data is not None:
                    return (deepcopy(entry.data), library, entry)
        return None

    def getSoluteDataFromGroups(self, species):
        """
        Return the set of Abraham solute parameters corresponding to a given
        :class:`Species` object `species` by estimation using the Platts group
        additivity method. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        
        It averages (linearly) over the desciptors for each Molecule (resonance isomer)
        in the Species.
        """       
        soluteData = SoluteData(0.0,0.0,0.0,0.0,0.0)
        count = 0
        comments = []
        for molecule in species.molecule:
            molecule.clearLabeledAtoms()
            molecule.updateAtomTypes()
            sdata = self.estimateSoluteViaGroupAdditivity(molecule)
            soluteData.S += sdata.S
            soluteData.B += sdata.B
            soluteData.E += sdata.E
            soluteData.L += sdata.L
            soluteData.A += sdata.A
            count += 1
            comments.append(sdata.comment)
        
        soluteData.S /= count
        soluteData.B /= count
        soluteData.E /= count
        soluteData.L /= count
        soluteData.A /= count
        
        # Print groups that are used for debugging purposes
        soluteData.comment = "Average of {0}".format(" and ".join(comments))

        return soluteData
   
    def transformLonePairs(self, molecule):
        """
        Changes lone pairs in a molecule to two radicals for purposes of finding
        solute data via group additivity. Transformed for each atom based on valency.
        """
        saturatedStruct = molecule.copy(deep=True)
        addedToPairs = {}

        for atom in saturatedStruct.atoms:
            addedToPairs[atom] = 0
            if atom.lonePairs > 0:
                charge = atom.charge # Record this so we can conserve it when checking
                bonds = saturatedStruct.getBonds(atom)
                sumBondOrders = 0
                for key, bond in bonds.iteritems():
                    if bond.order == 'S': sumBondOrders += 1
                    if bond.order == 'D': sumBondOrders += 2
                    if bond.order == 'T': sumBondOrders += 3
                    if bond.order == 'B': sumBondOrders += 1.5 # We should always have 2 'B' bonds (but what about Cbf?)
                if atomTypes['Val4'] in atom.atomType.generic: # Carbon, Silicon
                    while(atom.radicalElectrons + charge + sumBondOrders < 4):
                        atom.decrementLonePairs()
                        atom.incrementRadical()
                        atom.incrementRadical()
                        addedToPairs[atom] += 1
                if atomTypes['Val5'] in atom.atomType.generic: # Nitrogen
                    while(atom.radicalElectrons + charge + sumBondOrders < 3):
                        atom.decrementLonePairs()
                        atom.incrementRadical()
                        atom.incrementRadical()
                        addedToPairs[atom] += 1
                if atomTypes['Val6'] in atom.atomType.generic: # Oxygen, sulfur
                    while(atom.radicalElectrons + charge + sumBondOrders < 2):
                        atom.decrementLonePairs()
                        atom.incrementRadical()
                        atom.incrementRadical()
                        addedToPairs[atom] += 1
                if atomTypes['Val7'] in atom.atomType.generic: # Chlorine
                    while(atom.radicalElectrons + charge + sumBondOrders < 1):
                        atom.decrementLonePairs()
                        atom.incrementRadical()
                        atom.incrementRadical()
                        addedToPairs[atom] += 1

        saturatedStruct.update()
        saturatedStruct.updateLonePairs()
        
        return saturatedStruct, addedToPairs

    def removeHBonding(self, saturatedStruct, addedToRadicals, addedToPairs, soluteData):  
        
        # Remove hydrogen bonds and restore the radical
        for atom in addedToRadicals:
            for H, bond in addedToRadicals[atom]:
                saturatedStruct.removeBond(bond)
                saturatedStruct.removeAtom(H)
                atom.incrementRadical()

        # Change transformed lone pairs back
        for atom in addedToPairs:    
            if addedToPairs[atom] > 0:
                for pair in range(1, addedToPairs[atom]):
                    saturatedStruct.decrementRadical()
                    saturatedStruct.decrementRadical()
                    saturatedStruct.incrementLonePairs()

        # Update Abraham 'A' H-bonding parameter for unsaturated struct
        for atom in saturatedStruct.atoms:
            # Iterate over heavy (non-hydrogen) atoms
            if atom.isNonHydrogen() and atom.radicalElectrons > 0:
                for electron in range(1, atom.radicalElectrons):
                    # Get solute data for radical group    
                    try:
                        self.__addGroupSoluteData(soluteData, self.groups['radical'], saturatedStruct, {'*':atom})
                    except KeyError: pass
      
        return soluteData

    def estimateSoluteViaGroupAdditivity(self, molecule):
        """
        Return the set of Abraham solute parameters corresponding to a given
        :class:`Molecule` object `molecule` by estimation using the Platts' group
        additivity method. If no group additivity values are loaded, a
        :class:`DatabaseError` is raised.
        """
        # For thermo estimation we need the atoms to already be sorted because we
        # iterate over them; if the order changes during the iteration then we
        # will probably not visit the right atoms, and so will get the thermo wrong
        molecule.sortAtoms()

        # Create the SoluteData object with the intercepts from the Platts groups
        soluteData = SoluteData(
            S = 0.277,
            B = 0.071,
            E = 0.248,
            L = 0.13,
            A = 0.003
        )
        
        addedToRadicals = {} # Dictionary of key = atom, value = dictionary of {H atom: bond}
        addedToPairs = {} # Dictionary of key = atom, value = # lone pairs changed
        saturatedStruct = molecule.copy(deep=True)

        # Convert lone pairs to radicals, then saturate with H.
       
        # Change lone pairs to radicals based on valency
        if sum([atom.lonePairs for atom in saturatedStruct.atoms]) > 0: # molecule contains lone pairs
            saturatedStruct, addedToPairs = self.transformLonePairs(saturatedStruct)

        # Now saturate radicals with H
        if sum([atom.radicalElectrons for atom in saturatedStruct.atoms]) > 0: # radical species
            addedToRadicals = saturatedStruct.saturate()

        # Saturated structure should now have no unpaired electrons, and only "expected" lone pairs
        # based on the valency
        for atom in saturatedStruct.atoms:
            # Iterate over heavy (non-hydrogen) atoms
            if atom.isNonHydrogen():
                # Get initial solute data from main group database. Every atom must
                # be found in the main abraham database
                try:
                    self.__addGroupSoluteData(soluteData, self.groups['abraham'], saturatedStruct, {'*':atom})
                except KeyError:
                    logging.error("Couldn't find in main abraham database:")
                    logging.error(saturatedStruct)
                    logging.error(saturatedStruct.toAdjacencyList())
                    raise
                # Get solute data for non-atom centered groups (being found in this group
                # database is optional)    
                try:
                    self.__addGroupSoluteData(soluteData, self.groups['nonacentered'], saturatedStruct, {'*':atom})
                except KeyError: pass
        
        soluteData = self.removeHBonding(saturatedStruct, addedToRadicals, addedToPairs, soluteData)

        return soluteData

    def __addGroupSoluteData(self, soluteData, database, molecule, atom):
        """
        Determine the Platts group additivity solute data for the atom `atom`
        in the structure `structure`, and add it to the existing solute data
        `soluteData`.
        """

        node0 = database.descendTree(molecule, atom, None)

        if node0 is None:
            raise KeyError('Node not found in database.')

        # It's possible (and allowed) that items in the tree may not be in the
        # library, in which case we need to fall up the tree until we find an
        # ancestor that has an entry in the library
        node = node0
        
        while node is not None and node.data is None:
            node = node.parent
        if node is None:
            raise KeyError('Node has no parent with data in database.')
        data = node.data
        comment = node.label
        while isinstance(data, basestring) and data is not None:
            for entry in database.entries.values():
                if entry.label == data:
                    data = entry.data
                    comment = entry.label
                    break
        comment = '{0}({1})'.format(database.label, comment)

        # This code prints the hierarchy of the found node; useful for debugging
        #result = ''
        #while node is not None:
        #   result = ' -> ' + node + result
        #   node = database.tree.parent[node]
        #print result[4:]
        
        # Add solute data for each atom to the overall solute data for the molecule.
        soluteData.S += data.S
        soluteData.B += data.B
        soluteData.E += data.E
        soluteData.L += data.L
        soluteData.A += data.A
        soluteData.comment += comment + "+"
        
        return soluteData

    
    def calcH(self, soluteData, solventData):
        """
        Returns the enthalpy of solvation, at 298K, in J/mol
        """
        # Use Mintz parameters for solvents. Multiply by 1000 to go from kJ->J to maintain consistency
        delH = 1000*((soluteData.S*solventData.s_h)+(soluteData.B*solventData.b_h)+(soluteData.E*solventData.e_h)+(soluteData.L*solventData.l_h)+(soluteData.A*solventData.a_h)+solventData.c_h)  
        return delH
    
    def calcG(self, soluteData, solventData):
        """
        Returns the Gibbs free energy of solvation, at 298K, in J/mol
        """
        # Use Abraham parameters for solvents to get log K
        logK = (soluteData.S*solventData.s_g)+(soluteData.B*solventData.b_g)+(soluteData.E*solventData.e_g)+(soluteData.L*solventData.l_g)+(soluteData.A*solventData.a_g)+solventData.c_g
        # Convert to delG with units of J/mol
        delG = -8.314*298*2.303*logK
        return delG
        
    def calcS(self, delG, delH):
        """
        Returns the entropy of solvation, at 298K, in J/mol/K
        """
        delS = (delH-delG)/298
        return delS
    
    def getSolvationCorrection(self, soluteData, solventData):
        """ 
        Given a soluteData and solventData object, calculates the enthalpy, entropy,
        and Gibbs free energy of solvation at 298 K. Returns a SolvationCorrection
        object
        """
        correction = SolvationCorrection(0.0, 0.0, 0.0)
        correction.enthalpy = self.calcH(soluteData, solventData)
        correction.gibbs = self.calcG(soluteData, solventData)  
        correction.entropy = self.calcS(correction.gibbs, correction.enthalpy) 
        return correction
    
    def getSolvationKinetics(self, solvationReaction):
        """
        Return the barrier correction for a reaction
        """
        raise NotImplementedError()
