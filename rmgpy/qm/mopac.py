import os
import re
import external.cclib as cclib
import logging
from subprocess import Popen
import distutils.spawn

from rmgpy.molecule import Molecule
from qmdata import CCLibData
from molecule import QMMolecule


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
    
    failureKeys = ['IMAGINARY FREQUENCIES', 'EXCESS NUMBER OF OPTIMIZATION CYCLES', 'NOT ENOUGH TIME FOR ANOTHER CYCLE']
    
    def writeInputFile(self, attempt, top_keys, bottom_keys, polar_keys):
        """
        Check's that an output file exists and was successful.
        
        Returns a boolean flag that states whether a successful MOPAC simulation already exists for the molecule with the 
        given (augmented) InChI Key.
        """
        
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "mop")
        mol = openbabel.OBMol()
    
        if attempt <= scriptAttempts: #use UFF-refined coordinates
            obConversion.ReadFile(mol, self.geometry.getRefinedMolFilePath() )
        else:
            obConversion.ReadFile(mol, self.geometry.getCrudeMolFilePath() )
    
        mol.SetTitle(self.geometry.uniqueID) 
        obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
    
        input_string = obConversion.WriteString(mol)
        
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
    
    def getParser(self, outputFile):
        """
        Returns the appropriate cclib parser.
        """
        return cclib.parser.Mopac(outputFile)

class MopacMol(QMMolecule, Mopac):
    """
    A base Class for calculations of molecules using MOPAC. 
    
    Inherits from both :class:`QMMolecule` and :class:`Mopac`.
    """

    #: Keywords that will be added at the top and bottom of the qm input file
    keywords = [
                {'top':"precise nosym", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym gnorm=0.0 nonr", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym gnorm=0.0", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym gnorm=0.0 bfgs", 'bottom':"oldgeo thermo nosym precise "},
                {'top':"precise nosym recalc=10 dmax=0.10 nonr cycles=2000 t=2000", 'bottom':"oldgeo thermo nosym precise "},
                ]

    def writeInputFile(self, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attempt`.
        """
        
        molfile = self.getMolFilePathForCalculation(attempt) 
        atomline = re.compile('\s*([\- ][0-9.]+)\s+([\- ][0-9.]+)+\s+([\- ][0-9.]+)\s+([A-Za-z]+)')
        
        output = [ self.geometry.uniqueIDlong, '' ]
 
        atomCount = 0
        with open(molfile) as molinput:
            for line in molinput:
                match = atomline.match(line)
                if match:
                    output.append("{0:4s} {1} 1 {2} 1 {3} 1".format(match.group(4), match.group(1), match.group(2), match.group(3)))
                    atomCount += 1
        assert atomCount == len(self.molecule.atoms)
    
        output.append('')
        input_string = '\n'.join(output)
        
        top_keys, bottom_keys, polar_keys = self.inputFileKeywords(attempt)
        with open(self.inputFilePath, 'w') as mopacFile:
            mopacFile.write(top_keys)
            mopacFile.write('\n')
            mopacFile.write(input_string)
            mopacFile.write('\n')
            mopacFile.write(bottom_keys)
            if self.usePolar:
                mopacFile.write('\n\n\n')
                mopacFile.write(polar_keys)
                
    def inputFileKeywords(self, attempt):
        """
        Return the top, bottom, and polar keywords.
        """
        raise NotImplementedError("Should be defined by subclass, eg. MopacMolPM3")
        
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

class MopacMolPM3(MopacMolPMn):
    """
    Mopac PM3 calculations for molecules
    
    This is a class of its own in case you wish to do anything differently,
    but for now it's the same as all the MOPAC PMn calculations, only pm3
    """
    pm_method = 'pm3'

class MopacMolPM6(MopacMolPMn):
    """
    Mopac PM6 calculations for molecules
    
    This is a class of its own in case you wish to do anything differently,
    but for now it's the same as all the MOPAC PMn calculations, only pm6
    """
    pm_method = 'pm6'

class MopacMolPM7(MopacMolPMn):
    """
    Mopac PM7 calculations for molecules
    
    This is a class of its own in case you wish to do anything differently,
    but for now it's the same as all the MOPAC PMn calculations, only pm7
    """
    pm_method = 'pm7'
