import os
import re
import external.cclib as cclib
import logging
import time
from subprocess import Popen, PIPE
import distutils.spawn

from rmgpy.molecule import Molecule# Group
from rmgpy.reaction import Reaction
from qmdata import CCLibData
from molecule import QMMolecule
from reaction import QMReaction
from rmgpy.data.base import Entry
from rmgpy.data.kinetics import saveEntry
from rmgpy.data.kinetics.transitionstates import DistanceData


class Mopac:
    """
    A base class for all QM calculations that use MOPAC.
    
    Classes such as :class:`MopacMol` will inherit from this class.
    """

    inputFileExtension = '.mop'
    outputFileExtension = '.out'
    
    try:
        executablePath = distutils.spawn.find_executable('mopac') or \
                             distutils.spawn.find_executable('MOPAC2009.exe') or \
                             distutils.spawn.find_executable('MOPAC2012.exe')
    except:
        logging.debug("Did not find MOPAC on path, checking if it exists in a declared MOPAC_DIR...")
        mopacEnv = os.getenv('MOPAC_DIR', default="/opt/mopac")
        if os.path.exists(os.path.join(mopacEnv , 'MOPAC2012.exe')):
            executablePath = os.path.join(mopacEnv , 'MOPAC2012.exe')
        elif os.path.exists(os.path.join(mopacEnv , 'MOPAC2009.exe')):
            executablePath = os.path.join(mopacEnv , 'MOPAC2009.exe')
        elif os.path.exists(os.path.join(mopacEnv , 'mopac')):
            executablePath = os.path.join(mopacEnv , 'mopac')
        else:      
            executablePath = os.path.join(mopacEnv , '(MOPAC 2009 or 2012)')

    usePolar = False #use polar keyword in MOPAC
    
    "Keywords for the multiplicity"
    multiplicityKeywords = {}
    multiplicityKeywords[1] = ''
    multiplicityKeywords[2] = 'uhf doublet'
    multiplicityKeywords[3] = 'uhf triplet'
    multiplicityKeywords[4] = 'uhf quartet'
    multiplicityKeywords[5] = 'uhf quintet'
    multiplicityKeywords[6] = 'uhf sextet'
    multiplicityKeywords[7] = 'uhf septet'
    multiplicityKeywords[8] = 'uhf octet'
    multiplicityKeywords[9] = 'uhf nonet'
    
    "Keywords that will be added at the top of the qm input file"
    keywordsTop = {}
    keywordsTop[1] = "precise nosym"
    keywordsTop[2] = "precise nosym gnorm=0.0 nonr"
    keywordsTop[3] = "precise nosym gnorm=0.0"
    keywordsTop[4] = "precise nosym gnorm=0.0 bfgs"
    keywordsTop[5] = "precise nosym recalc=10 dmax=0.10 nonr cycles=2000 t=2000"
    
    "Keywords that will be added at the bottom of the qm input file"
    keywordsBottom = {}
    keywordsBottom[1] = "oldgeo thermo nosym precise "
    keywordsBottom[2] = "oldgeo thermo nosym precise "
    keywordsBottom[3] = "oldgeo thermo nosym precise "
    keywordsBottom[4] = "oldgeo thermo nosym precise "
    keywordsBottom[5] = "oldgeo thermo nosym precise "
    
    scriptAttempts = len(keywordsTop)
    maxAttempts = 2 * scriptAttempts
    
        return self.verifyOutputFile()
        
    def parse(self):
        """
        Parses the results of the Mopac calculation, and returns a CCLibData object.
        """
        parser = cclib.parser.Mopac(self.outputFilePath)
        parser.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
        cclibData = parser.parse()
        radicalNumber = sum([i.radicalElectrons for i in self.molecule.atoms])
        qmData = CCLibData(cclibData, radicalNumber+1)
        return qmData

class MopacMol(QMMolecule, Mopac):
    """
    A base Class for calculations of molecules using MOPAC. 
    
    Inherits from both :class:`QMMolecule` and :class:`Mopac`.
    """
                
    def inputFileKeywords(self, attempt):
        """
        Return the top, bottom, and polar keywords.
        """
        raise NotImplementedError("Should be defined by subclass, eg. MopacMolPM3")
        
    def writeInputFile(self, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "mop")
        mol = openbabel.OBMol()
    
        obConversion.ReadFile(mol, self.getMolFilePathForCalculation(attempt) )
        
        mol.SetTitle(self.geometry.uniqueIDlong)
        obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
        input_string = obConversion.WriteString(mol)
        top_keys, bottom_keys, polar_keys = self.inputFileKeywords(attempt)
        with open(self.inputFilePath, 'w') as mopacFile:
            mopacFile.write(top_keys)
            mopacFile.write(input_string)
            mopacFile.write('\n')
            mopacFile.write(bottom_keys)
            if self.usePolar:
                mopacFile.write('\n\n\n')
                mopacFile.write(polar_keys)
    
    def verifyOutputFile(self):
        """
        Check's that an output file exists and was successful.
        
        Returns a boolean flag that states whether a successful MOPAC simulation already exists for the molecule with the 
        given (augmented) InChI Key.
        """
        
        if not os.path.exists(self.outputFilePath):
            logging.info("Output file {0} does not exist.".format(self.outputFilePath))
            return False
        
        InChIMatch=False #flag (1 or 0) indicating whether the InChI in the file matches InChIaug this can only be 1 if InChIFound is also 1
        InChIFound=False #flag (1 or 0) indicating whether an InChI was found in the log file
        
        with open(inputFilePath, 'w') as mopacFile:
            mopacFile.write(top_keys)
            mopacFile.write(input_string)
            mopacFile.write('\n')
            mopacFile.write(bottom_keys)
            if self.usePolar:
                mopacFile.write('\n\n\n')
                mopacFile.write(polar_keys)
        
        return self.geometry.uniqueID + self.inputFileExtension
    
    def convertMolFile(self, outputType, attempt, scriptAttempts):
        if attempt <= scriptAttempts: #use UFF-refined coordinates
            inputFilePath = self.geometry.getRefinedMolFilePath()
        else:
            inputFilePath = self.geometry.getCrudeMolFilePath()
                    
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", outputType)
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, inputFilePath)
        mol.SetTitle(self.geometry.uniqueID) 
        obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
        input_string = obConversion.WriteString(mol)
        
        return input_string
       
    def run(self, inputFileName):
        # submits the input file to mopac
        command = os.path.join(self.directory, self.geometry.uniqueID + self.inputFileExtension)
        process = Popen([self.executablePath, command])
        process.communicate()# necessary to wait for executable termination!
    
        return self.checkNoFailure()
        
    def checkNoFailure(self):
        """
        checks whether the output file contains any of the 
        failure keywords
        """
        file = os.path.join(self.directory,self.geometry.uniqueID+self.outputFileExtension)
        with open(file) as qmfile:    
            for each_line in qmfile:
                each_line = each_line.rstrip().strip()
                for element in self.failureKeys:#search for failure keywords
                    if element in each_line:
                        logging.error("MOPAC output file contains the following error %s")%element
                        return False
                
                for element in self.successKeys: #search for success keywords
                    if element in line:
                        successKeysFound[element] = True
                
                if "InChI=" in line:
                    logFileInChI = line #output files should take up to 240 characters of the name in the input file
                    InChIFound = True
                    if self.uniqueIDlong in logFileInChI:
                        InChIMatch = True
                    elif self.uniqueIDlong.startswith(logFileInChI):
                        logging.info("InChI too long to check, but beginning matches so assuming OK.")
                        InChIMatch = True
                    else:
                        logging.warning("InChI in log file ({0}) didn't match that in geometry ({1}).".format(logFileInChI, self.uniqueIDlong))                    
                        # Use only up to first 80 characters to match due to MOPAC bug which deletes 81st character of InChI string
                        if self.uniqueIDlong.startswith(logFileInChI[:80]):
                            logging.warning("but the beginning matches so it's probably just a truncation problem.")
                            InChIMatch = True
        # Check that ALL 'success' keywords were found in the file.
        if not all( successKeysFound.values() ):
            logging.error('Not all of the required keywords for success were found in the output file!')
            return False
        
        if not InChIFound:
            logging.error("No InChI was found in the MOPAC output file {0}".format(self.outputFilePath))
            return False
        
        if not InChIMatch:
            #InChIs do not match (most likely due to limited name length mirrored in log file (240 characters), but possibly due to a collision)
            return self.checkForInChiKeyCollision(logFileInChI) # Not yet implemented!

        # Compare the optimized geometry to the original molecule
        qmData = self.parse()
        cclibMol = Molecule()
        cclibMol.fromXYZ(qmData.atomicNumbers, qmData.atomCoords.value)
        testMol = self.molecule.toSingleBonds()
        if not cclibMol.isIsomorphic(testMol):
            logging.info("Incorrect connectivity for optimized geometry in file {0}".format(self.outputFilePath))
            return False

        logging.info("Successful {1} quantum result in {0}".format(self.outputFilePath, self.__class__.__name__))
        return True
        
        #InChIs do not match (most likely due to limited name length mirrored in log file (240 characters), but possibly due to a collision)
        return self.checkForInChiKeyCollision(logFileInChI) # Not yet implemented!
            
    def generateQMData(self):
        """
        Calculate the QM data and return a QMData object, or None if it fails.
        """
        for atom in self.molecule.vertices:
            if atom.atomType.label in ('N5s', 'N5d', 'N5dd', 'N5t', 'N5b'):
                return None

        if self.verifyOutputFile():
            logging.info("Found a successful output file already; using that.")
            source = "QM {0} calculation found from previous run.".format(self.__class__.__name__)
        else:
            self.createGeometry()
            success = False
            for attempt in range(1, self.maxAttempts+1):
                self.writeInputFile(attempt)
                logging.info('Trying {3} attempt {0} of {1} on molecule {2}.'.format(attempt, self.maxAttempts, self.molecule.toSMILES(), self.__class__.__name__))
                success = self.run()
                if success:
                    logging.info('Attempt {0} of {1} on species {2} succeeded.'.format(attempt, self.maxAttempts, self.molecule.toAugmentedInChI()))
                    source = "QM {0} calculation attempt {1}".format(self.__class__.__name__, attempt )
                    break
            else:
                logging.error('QM thermo calculation failed for {0}.'.format(self.molecule.toAugmentedInChI()))
                return None
        result = self.parse() # parsed in cclib
        result.source = source
        return result # a CCLibData object


class MopacMolPMn(MopacMol):
    """
    Mopac PMn calculations for molecules (n undefined here)
    
    This is a parent class for MOPAC PMn calculations.
    Inherit it, and define the pm_method, then redefine 
    anything you wish to do differently.
    """
    pm_method = '(should be defined by sub class)'
    def inputFileKeywords(self, attempt):
        """
        Return the top, bottom, and polar keywords for attempt number `attempt`.
        
        NB. `attempt` begins at 1, not 0.
        """
        assert attempt <= self.maxAttempts
        
        if attempt > self.scriptAttempts:
            attempt -= self.scriptAttempts
        
        multiplicity_keys = self.multiplicityKeywords[self.geometry.molecule.multiplicity]

        top_keys = "{method} {mult} {top}".format(
                method = self.pm_method,
                mult = multiplicity_keys,
                top = self.keywords[attempt-1]['top'],
                )
        bottom_keys = "{bottom} {method} {mult}".format(
                method = self.pm_method,
                bottom = self.keywords[attempt-1]['bottom'],
                mult = multiplicity_keys,
                )
        polar_keys = "oldgeo {polar} nosym precise {method} {mult}".format(
                method = self.pm_method,
                polar = ('polar' if self.geometry.molecule.multiplicity == 1 else 'static'),
                mult = multiplicity_keys,
                )

        return top_keys, bottom_keys, polar_keys

#TS
class MopacTS(QMReaction, Mopac):
    #*****change this for TS
    "Keywords for the multiplicity"
    multiplicityKeywords = {}
    multiplicityKeywords[1] = ''
    multiplicityKeywords[2] = 'uhf doublet'
    multiplicityKeywords[3] = 'uhf triplet'
    multiplicityKeywords[4] = 'uhf quartet'
    multiplicityKeywords[5] = 'uhf quintet'
    multiplicityKeywords[6] = 'uhf sextet'
    multiplicityKeywords[7] = 'uhf septet'
    multiplicityKeywords[8] = 'uhf octet'
    multiplicityKeywords[9] = 'uhf nonet'

    "Keywords that will be added at the top of the qm input file"
    keywordsTop = {}
    keywordsTop[1] = "ts"
    keywordsTop[2] = "ts recalc=5"
    keywordsTop[3] = "ts ddmin=0.0001"
    keywordsTop[4] = "ts recalc=5 ddmin=0.0001"

    "Keywords that will be added at the bottom of the qm input file"
    keywordsBottom = {}
    keywordsBottom[1] = "oldgeo force"
    keywordsBottom[2] = "oldgeo force esp"
    keywordsBottom[3] = "oldgeo force vectors"
    keywordsBottom[4] = "oldgeo force vectors esp"

    scriptAttempts = len(keywordsTop)

    failureKeys = ['GRADIENT IS TOO LARGE', 
                'EXCESS NUMBER OF OPTIMIZATION CYCLES', 
                'NOT ENOUGH TIME FOR ANOTHER CYCLE',
                '6 IMAGINARY FREQUENCIES',
                '5 IMAGINARY FREQUENCIES',
                '4 IMAGINARY FREQUENCIES',
                '3 IMAGINARY FREQUENCIES',
                '2 IMAGINARY FREQUENCIES'
                ]
    
    def runIRC(self):
        self.testReady()
        # submits the input file to mopac
        process = Popen([self.executablePath, self.inputFilePath])
        process.communicate()# necessary to wait for executable termination!
    
        return self.verifyIRCOutputFile()
    
    def writeInputFile(self, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "mop")
        mol = openbabel.OBMol()
    
        obConversion.ReadFile(mol, self.geometry.getRefinedMolFilePath() )
        
        mol.SetTitle(self.uniqueID)
        obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
        input_string = obConversion.WriteString(mol)
        top_keys, bottom_keys, polar_keys = self.inputFileKeywords(attempt, 2)
        with open(self.inputFilePath, 'w') as mopacFile:
            mopacFile.write(top_keys)
            mopacFile.write(input_string)
            mopacFile.write('\n')
            mopacFile.write(bottom_keys)
            if self.usePolar:
                mopacFile.write('\n\n\n')
                mopacFile.write(polar_keys)
                
    def writeIRCFile(self):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "mop")
        parseOutput = cclib.parser.Mopac(self.outputFilePath.split('IRC')[0] + '.out')
        parseOutput = parseOutput.parse()
        reload(openbabel)
        mol = cclib.bridge.makeopenbabel(parseOutput.atomcoords[0], parseOutput.atomnos)
        mol.SetTitle(self.geometry.uniqueIDlong)
        obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
        input_string = obConversion.WriteString(mol)
    
        top_keys = 'irc=1*'
        with open(self.inputFilePath, 'w') as mopacFile:
            mopacFile.write(top_keys)
            mopacFile.write(input_string)
            mopacFile.write('\n')
    
    def verifyOutputFile(self):
        """
        Check's that an output file exists and was successful.
        
        Returns a boolean flag that states whether a successful MOPAC simulation already exists for the molecule with the 
        given (augmented) InChI Key.
        """
        
        if not os.path.exists(self.outputFilePath):
            logging.info("Output file {0} does not exist.".format(self.outputFilePath))
            return False, False
        
        # Initialize dictionary with "False"s 
        successKeysFound = dict([(key, False) for key in self.successKeys])
        
        with open(self.outputFilePath) as outputFile:
            for line in outputFile:
                line = line.strip()
                
                for element in self.failureKeys: #search for failure keywords
                    if element in line:
                        logging.error("MOPAC output file contains the following error: {0}".format(element) )
                        return False, False
                    
                for element in self.successKeys: #search for success keywords
                    if element in line:
                        successKeysFound[element] = True
               
        # Check that ALL 'success' keywords were found in the file.
        if not successKeysFound['MOPAC DONE']:
            logging.error('Not all of the required keywords for success were found in the output file!')
            return False, False
        else:
            logging.info("Successful MOPAC quantum result found in {0}".format(self.outputFilePath))
            return True, False
        
        #InChIs do not match (most likely due to limited name length mirrored in log file (240 characters), but possibly due to a collision)
        return self.checkForInChiKeyCollision(logFileInChI) # Not yet implemented!
    
    def convertMol(self, geomLines):
        atomcoords = []
        atomnos = []
        for line in geomLines:
            atType, x, y, z = line.split()
            if atType == 'H':
                atNum = 1
            elif atType == 'C':
                atNum = 6
            elif atType == 'O':
                atNum = 8
            coords = [float(x),float(y),float(z)]
            atomnos.append(atNum)
            atomcoords.append(coords)
        atomnos = numpy.array(atomnos, dtype=int)
        atomcoords = numpy.array(atomcoords)
        reload(openbabel)
        mol = cclib.bridge.makeopenbabel(atomcoords, atomnos) 
    
        return Molecule().fromOBMol(mol)
    
    def verifyIRCOutputFile(self):
        """
        Check's that an output file exists and was successful.
        
        Returns a boolean flag that states whether a successful MOPAC simulation already exists for the molecule with the 
        given (augmented) InChI Key.
        """
        
        if not os.path.exists(self.outputFilePath):
            logging.info("Output file {0} does not exist.".format(self.outputFilePath))
            return False
        
        # Initialize dictionary with "False"s 
        successKeysFound = dict([(key, False) for key in self.successKeys])
        
        with open(self.outputFilePath) as outputFile:
            for line in outputFile:
                line = line.strip()
                
                for element in self.failureKeys: #search for failure keywords
                    if element in line:
                        logging.error("MOPAC output file contains the following error: {0}".format(element) )
                        return False
                    
                for element in self.successKeys: #search for success keywords
                    if element in line:
                        successKeysFound[element] = True
        
        if not successKeysFound['MOPAC DONE']:
            logging.error('Not all of the required keywords for success were found in the IRC output file!')
            return False
        
        with open(self.outputFilePath.split('.')[0] + '.xyz') as geomFile:
            geomFile.pop(0)
            geomFile.pop(0)
            geom1 = []
            for line in geomFile:
                if not line.startswith('  reversed'):
                    geom1.append(line)
                else:
                    break
            geom1.pop()
            
            geom2 = []
            for line in reversed(geomFile):
                if not line.startswith(' DRC'):
                    geom2.append(line)
                else:
                    break
            geom2.pop()
        
        mol1 = self.convertMol(geom1)
        mol2 = self.convertMol(geom2)
        
        targetReaction = Reaction(
                                reactants = [reactant.toSingleBonds() for reactant in self.reaction.reactants],
                                products = [product.toSingleBonds() for product in self.reaction.products],
                                )
        testReaction = Reaction(
                                reactants = mol1.split(),
                                products = mol2.split(),                     
                                )
                                
        if targetReaction.isIsomorphic(testReaction):
            return True
        else:
            return False
    
    def parseTS(self, labels):
    
        tsParse = cclib.parser.Mopac(os.path.join(self.file_store_path, self.uniqueID + '.log'))
        tsParse = tsParse.parse()
    
        atom1 = openbabel.OBAtom()
        atom2 = openbabel.OBAtom()
        atom3 = openbabel.OBAtom()
    
        atom1.SetAtomicNum(int(tsParse.atomnos[labels[0]]))
        atom2.SetAtomicNum(int(tsParse.atomnos[labels[1]]))
        atom3.SetAtomicNum(int(tsParse.atomnos[labels[2]]))
    
        atom1coords = tsParse.atomcoords[-1][labels[0]].tolist()
        atom2coords = tsParse.atomcoords[-1][labels[1]].tolist()
        atom3coords = tsParse.atomcoords[-1][labels[2]].tolist()
    
        atom1.SetVector(*atom1coords)
        atom2.SetVector(*atom2coords)
        atom3.SetVector(*atom3coords)
        
        # from rmgpy.molecule.element import getElement
        # at1 = getElement(atom1.GetAtomicNum()).symbol
        # at2 = getElement(atom2.GetAtomicNum()).symbol
        # at3 = getElement(atom3.GetAtomicNum()).symbol
    
        atomDist = [str(atom1.GetDistance(atom2)), str(atom2.GetDistance(atom3)), str(atom1.GetDistance(atom3))]
    
        return atomDist
    
    def writeRxnOutputFile(self, labels):
        
        product = self.reaction.products[0].merge(self.reaction.products[1])
        star3 = product.getLabeledAtom('*1').sortingLabel
        star1 = product.getLabeledAtom('*3').sortingLabel
        product.atoms[star1].label = '*1'
        product.atoms[star3].label = '*3'
        
        atomDist = self.parseTS(labels)
        
        distances = {'d12':float(atomDist[0]), 'd23':float(atomDist[1]), 'd13':float(atomDist[2])}
        user = "Pierre Bhoorasingh <bhoorasingh.p@husky.neu.edu>"
        description = "Found via group estimation strategy using automatic transition state generator"
        entry = Entry(
            index = 1,
            item = self.reaction,
            data = DistanceData(distances=distances, method='B3LYP/6-31+G(d,p)'),
            shortDesc = "B3LYP/6-31+G(d,p) calculation via group estimated TS generator.",
            history = [(time.asctime(), user, 'action', description)]
        )
        
        outputDataFile = os.path.join(self.file_store_path, self.uniqueID + '.data')
        with open(outputDataFile, 'w') as parseFile:
            saveEntry(parseFile, entry)
    
class MopacTSPM3(MopacTS):
    def inputFileKeywords(self, attempt, multiplicity):
        """
        Inherits the writeInputFile methods from mopac.py
        """
        multiplicity_keys = self.multiplicityKeywords[multiplicity]

        top_keys = "pm3 {0} {1}".format(
                multiplicity_keys,
                self.keywordsTop[attempt],
                )
        bottom_keys = "{0} pm3 {1}".format(
                self.keywordsBottom[attempt],
                multiplicity_keys,
                )
        polar_keys = "oldgeo {0} nosym precise pm3 {1}".format(
                'polar' if multiplicity == 1 else 'static',
                multiplicity_keys,
                )

        return top_keys, bottom_keys, polar_keys

class MopacTSPM7(MopacTS):
    def inputFileKeywords(self, attempt, multiplicity):
        """
        Inherits the writeInputFile methods from mopac.py
        """
        multiplicity_keys = self.multiplicityKeywords[multiplicity]

        top_keys = "pm7 {0} {1}".format(
                multiplicity_keys,
                self.keywordsTop[attempt],
                )
        bottom_keys = "{0} pm7 {1}".format(
                self.keywordsBottom[attempt],
                multiplicity_keys,
                )
        polar_keys = "oldgeo {0} nosym precise pm7 {1}".format(
                'polar' if multiplicity == 1 else 'static',
                multiplicity_keys,
                )

        return top_keys, bottom_keys, polar_keys
