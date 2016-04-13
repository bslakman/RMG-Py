#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from unittest import TestCase, TestLoader, TextTestRunner
from external.wip import work_in_progress

from rmgpy import settings
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.data.solvation import SolvationDatabase, SolvationKinetics
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.reaction import Reaction
from rmgpy.kinetics.diffusionLimited import LiquidKinetics, diffusionLimiter
###################################################

class TestDiffusionLimited(TestCase):
    
    def setUp(self):
        self.database = SolvationDatabase()
        self.database.load(os.path.join(settings['database.directory'], 'solvation'))
        self.barrierDatabase = SolvationKinetics()
        self.kineticsDatabase = KineticsDatabase()
        self.kineticsDatabase.load(os.path.join(settings['database.directory'], 'kinetics'), families=['H_Abstraction'], libraries=[])
        self.barrierDatabase.family = self.kineticsDatabase.families['H_Abstraction']
        self.barrierDatabase.load(os.path.join(settings['database.directory'], 'kinetics', 'families', 'H_Abstraction'), None, None)

        solventData = self.database.getSolventData('n-octane')
        diffusionLimiter.enable(solventData, self.database, self.kineticsDatabase)

    def runTest(self):
        pass
    
    def testCorrectIntrinsicRate(self):
        "Test we can correct a gas phase reaction rate"
        r1 = Molecule().fromAdjacencyList(
"""
1    C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {1,S}
"""
        ) # methane
        r2 = Molecule().fromAdjacencyList(
"""
1    O u1 p2 c0 {2,S}
2    H u0 p0 c0 {1,S}
"""
        ) # OH
        p1 = Molecule().fromAdjacencyList(
"""
1    C u1 p0 c0 {2,S} {3,S} {4,S}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
"""
        ) # methyl
        p2 = Molecule().fromAdjacencyList(
"""
1    O u0 p2 c0 {2,S} {3,S}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
"""
        ) # water
        reaction = Reaction(reactants=[r1,r2], products=[p1,p2]) 
        liquidKinetics = LiquidKinetics(reaction, diffusionLimiter.solventData, diffusionLimiter.kinetics)
        try:
            liquidKinetics.correctIntrinsicRate(500)
        except:
            self.fail("Failed to correct intrinsic rate")

#####################################################

if __name__ == '__main__':
    suite = TestLoader().loadTestsFromTestCase(TestDiffusionLimited)
    TextTestRunner(verbosity=2).run(suite)
