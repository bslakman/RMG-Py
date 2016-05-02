import unittest
import os.path
import shutil

import rmgpy
from rmgpy.rmg.main import RMG
from rmgpy.tools.generate_reactions import *

class GenerateReactionsTest(unittest.TestCase):

    def test(self):
        folder = os.path.join(os.path.dirname(rmgpy.__file__),'tools/data/generate')
        
        inputFile = os.path.join(folder,'input.py')
        
        rmg = RMG()
        rmg = execute(rmg, inputFile, folder)

        self.assertIsNotNone(rmg)
        self.assertIsNotNone(rmg.reactionModel.outputSpeciesList)
        self.assertIsNotNone(rmg.reactionModel.outputReactionList)


        shutil.rmtree(os.path.join(folder,'pdep'))

    def testDuplicateReaction(self):
        """
        Test that the radical addition reaction

        HCJ=O + CH2O = [CH2]OC=O

        present in the reaction library "Methylformate",
        only appears once in the model.

        """

        from rmgpy.reaction import Reaction
        from rmgpy.molecule import Molecule
        folder = os.path.join(os.path.dirname(rmgpy.__file__),'tools/data/generate/duplicates')
        
        inputFile = os.path.join(folder,'input.py')

        rmg = RMG()
        rmg = execute(rmg, inputFile, folder)

        self.assertIsNotNone(rmg)
        
        rxnFlagged = Reaction(reactants=[Molecule(SMILES='[CH]=O'),Molecule(SMILES='C=O')],
                       products=[Molecule(SMILES='[CH2]OC=O')])

        count = 0
        for reaction in rmg.reactionModel.core.reactions:
            if reaction.isIsomorphic(rxnFlagged):
                count += 1

        self.assertEquals(count, 1)

        shutil.rmtree(os.path.join(folder,'pdep'))
        
    def testLibraryReactionEntersCore(self):
        """
        Test that a reaction from a Reaction Library enters the core
        right after the initialization step if all the input species are 
        present in that reaction.
        
        The following reaction from the Methylformate library
        
        HCjO + CH2O <=> Fmoml
        
        should appear in the model if HCjO, CH2O and Fmoml are all used as input species
        """
        from rmgpy.reaction import Reaction
        from rmgpy.molecule import Molecule
        folder = os.path.join(os.path.dirname(rmgpy.__file__),'tools/data/generate/libraryReaction')
        
        inputFile = os.path.join(folder,'input.py')

        rmg = RMG()
        rmg = execute(rmg, inputFile, folder)

        self.assertIsNotNone(rmg)
        
        # Assert that the flagged reaction occurs
        rxnFlagged = Reaction(reactants=[Molecule(SMILES='[CH]=O'),Molecule(SMILES='C=O')],
                       products=[Molecule(SMILES='[CH2]OC=O')])

        count = 0
        for reaction in rmg.reactionModel.core.reactions:
            if reaction.isIsomorphic(rxnFlagged):
                count += 1

        self.assertEquals(count, 1)
        

        # Assert that the core only has 1 reaction
        self.assertEquals(len(rmg.reactionModel.core.reactions),1)
        shutil.rmtree(os.path.join(folder,'pdep'))
        
    
    def setUp(self):
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None
        
    def tearDown(self):
        """
        Reset the loaded database
        """
        import rmgpy.data.rmg
        rmgpy.data.rmg.database = None
