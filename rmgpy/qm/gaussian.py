import os

import external.cclib as cclib
import logging
import time
from subprocess import Popen, PIPE
import re

from rmgpy.molecule import Group, Molecule
from rmgpy.reaction import Reaction
from qmdata import CCLibData
from molecule import QMMolecule
from reaction import QMReaction
from rmgpy.data.base import Entry
from rmgpy.data.kinetics import saveEntry
from rmgpy.data.kinetics.transitionstates import DistanceData

class Gaussian:
    """
    A base class for all QM calculations that use Gaussian.
    
    Classes such as :class:`GaussianMol` will inherit from this class.
    """
    
    inputFileExtension = '.gjf'
    outputFileExtension = '.log'
    
    gaussEnv = os.getenv('GAUSS_EXEDIR') or os.getenv('g09root') or os.getenv('g03root') or ""
    
    # GAUSS_EXEDIR may be a list like "path1:path2:path3"
    for possibleDir in gaussEnv.split(':'):
        if os.path.exists(os.path.join(possibleDir , 'g09')):
            executablePath = os.path.join(possibleDir , 'g09')
            break
        elif os.path.exists(os.path.join(possibleDir , 'g03')):
            executablePath = os.path.join(possibleDir , 'g03')
            break
    else:
        executablePath = os.path.join(gaussEnv , '(g03 or g09)')

    usePolar = False
     
    def testReady(self):
        if not os.path.exists(self.executablePath):
            raise Exception("Couldn't find Gaussian executable at {0}. Try setting your GAUSS_EXEDIR environment variable.".format(self.executablePath))

   
    def run(self):
        self.testReady()
        # submits the input file to Gaussian
        process = Popen([self.executablePath, self.inputFilePath, self.outputFilePath])
        process.communicate()# necessary to wait for executable termination!
        
        return self.verifyOutputFile()
        
    def parse(self):
        """
        Parses the results of the Gaussian calculation, and returns a CCLibData object.
        """
        parser = cclib.parser.Gaussian(self.outputFilePath)
        parser.logger.setLevel(logging.ERROR) #cf. http://cclib.sourceforge.net/wiki/index.php/Using_cclib#Additional_information
        cclibData = parser.parse()
        radicalNumber = sum([i.radicalElectrons for i in self.molecule.atoms])
        qmData = CCLibData(cclibData, radicalNumber+1)
        return qmData
    
    
    
class GaussianMol(QMMolecule, Gaussian):
    """
    A base Class for calculations of molecules using Gaussian. 
    
    Inherits from both :class:`QMMolecule` and :class:`Gaussian`.
    """
    #: List of phrases to indicate success.
    #: ALL of these must be present in a successful job.
    successKeys = [
                   'Normal termination of Gaussian',
                  ]
    
    #: List of phrases that indicate failure
    #: NONE of these must be present in a succesful job.
    failureKeys = [
                   'ERROR TERMINATION',
                   'IMAGINARY FREQUENCIES'
                   ]
    
    def writeInputFile(self, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
    
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "gjf")
        mol = openbabel.OBMol()
    
        obConversion.ReadFile(mol, self.getMolFilePathForCalculation(attempt) )
    
        mol.SetTitle(self.geometry.uniqueIDlong)
        obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
        input_string = obConversion.WriteString(mol)
        top_keys = self.inputFileKeywords(attempt)
        with open(self.inputFilePath, 'w') as gaussianFile:
            gaussianFile.write(top_keys)
            gaussianFile.write(input_string)
            gaussianFile.write('\n')
            if self.usePolar:
                gaussianFile.write('\n\n\n')
                gaussianFile.write(polar_keys)
    
    def inputFileKeywords(self, attempt):
        """
        Return the top keywords.
        """
        raise NotImplementedError("Should be defined by subclass, eg. GaussianMolPM3")
    
    def generateQMData(self):
        """
        Calculate the QM data and return a QMData object.
        """
        self.createGeometry()
        if self.verifyOutputFile():
            logging.info("Found a successful output file already; using that.")
        else:
            success = False
            for attempt in range(1, self.maxAttempts+1):
                self.writeInputFile(attempt)
                success = self.run()
                if success:
                    logging.info('Attempt {0} of {1} on species {2} succeeded.'.format(attempt, self.maxAttempts, self.molecule.toAugmentedInChI()))
                    break
            else:
                raise Exception('QM thermo calculation failed for {0}.'.format(self.molecule.toAugmentedInChI()))
        result = self.parse() # parsed in cclib
        return result
    
    def verifyOutputFile(self):
        """
        Check's that an output file exists and was successful.
        
        Returns a boolean flag that states whether a successful GAUSSIAN simulation already exists for the molecule with the 
        given (augmented) InChI Key.
        
        The definition of finding a successful simulation is based on these criteria:
        1) finding an output file with the file name equal to the InChI Key
        2) NOT finding any of the keywords that are denote a calculation failure
        3) finding all the keywords that denote a calculation success.
        4) finding a match between the InChI of the given molecule and the InchI found in the calculation files
        
        If any of the above criteria is not matched, False will be returned and the procedures to start a new calculation 
        will be initiated.
        """
        if not os.path.exists(self.outputFilePath):
            logging.info("Output file {0} does not exist.".format(self.outputFilePath))
            return False
    
        InChIMatch=False #flag (1 or 0) indicating whether the InChI in the file matches InChIaug this can only be 1 if InChIFound is also 1
        InChIFound=False #flag (1 or 0) indicating whether an InChI was found in the log file
        
        # Initialize dictionary with "False"s 
        successKeysFound = dict([(key, False) for key in self.successKeys])
        
        with open(self.outputFilePath) as outputFile:
            for line in outputFile:
                line = line.strip()
                
                for element in self.failureKeys: #search for failure keywords
                    if element in line:
                        logging.error("Gaussian output file contains the following error: {0}".format(element) )
                        return False
                    
                for element in self.successKeys: #search for success keywords
                    if element in line:
                        successKeysFound[element] = True
               
                if line.startswith("InChI="):
                    logFileInChI = line #output files should take up to 240 characters of the name in the input file
                    InChIFound = True
                    if logFileInChI == self.geometry.uniqueIDlong:
                        InChIMatch = True
                    elif self.geometry.uniqueIDlong.startswith(logFileInChI):
                        logging.info("InChI too long to check, but beginning matches so assuming OK.")
                        InChIMatch = True
                    else:
                        logging.info("InChI in log file didn't match that in geometry.")
                        logging.info(self.geometry.uniqueIDlong)
                        logging.info(logFileInChI)
        
        # Check that ALL 'success' keywords were found in the file.
        if not all( successKeysFound.values() ):
            logging.error('Not all of the required keywords for success were found in the output file!')
            return False
        
        if not InChIFound:
            logging.error("No InChI was found in the Gaussian output file {0}".format(self.outputFilePath))
            return False
        
        if InChIMatch:
            logging.info("Successful Gaussian quantum result found in {0}".format(self.outputFilePath))
            # " + self.molfile.name + " ("+self.molfile.InChIAug+") has been found. This log file will be used.")
            return True
        else:
            return False # until the next line works
        
        #InChIs do not match (most likely due to limited name length mirrored in log file (240 characters), but possibly due to a collision)
        return self.checkForInChiKeyCollision(logFileInChI) # Not yet implemented!

class GaussianMolPM3(GaussianMol):

    #: Keywords that will be added at the top of the qm input files
    keywords = [
               "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3)",
               "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3) IOP(4/21=2)",
               "# pm3 opt=(verytight,calcfc,maxcyc=200) freq IOP(2/16=3) nosymm" ,
               "# pm3 opt=(verytight,calcfc,maxcyc=200) freq=numerical IOP(2/16=3) nosymm",
               "# pm3 opt=(verytight,gdiis,small) freq IOP(2/16=3)",
               "# pm3 opt=(verytight,nolinear,calcfc,small) freq IOP(2/16=3)",
               "# pm3 opt=(verytight,gdiis,maxcyc=200) freq=numerical IOP(2/16=3)",
               "# pm3 opt=tight freq IOP(2/16=3)",
               "# pm3 opt=tight freq=numerical IOP(2/16=3)",
               "# pm3 opt=(tight,nolinear,calcfc,small,maxcyc=200) freq IOP(2/16=3)",
               "# pm3 opt freq IOP(2/16=3)",
               "# pm3 opt=(verytight,gdiis) freq=numerical IOP(2/16=3) IOP(4/21=200)",
               "# pm3 opt=(calcfc,verytight,newton,notrustupdate,small,maxcyc=100,maxstep=100) freq=(numerical,step=10) IOP(2/16=3) nosymm",
               "# pm3 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",
               "# pm3 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",
               "# pm3 opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm",
               "# pm3 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm",
               "# pm3 opt=(calcall,small,maxcyc=100) IOP(2/16=3)",
               ]

    @property
    def scriptAttempts(self):
        "The number of attempts with different script keywords"
        return len(self.keywords)

    @property
    def maxAttempts(self):
        "The total number of attempts to try"
        return 2 * len(self.keywords)

    def inputFileKeywords(self, attempt):
        """
        Return the top keywords for attempt number `attempt`.

        NB. `attempt`s begin at 1, not 0.
        """
        assert attempt <= self.maxAttempts
        if attempt > self.scriptAttempts:
            attempt -= self.scriptAttempts
        return self.keywords[attempt-1]

class GaussianMolPM3(GaussianMol):
    """
    Gaussian PM3 calculations for molecules
    
    This is a class of its own in case you wish to do anything differently,
    but for now it's only the 'pm3' in the keywords that differs.
    """
    #: Keywords that will be added at the top of the qm input file
    keywords = [
                # The combinations of keywords were derived by Greg Magoon for pm3 in Gaussian. His comments are attached to each combination.
               "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3)", # added IOP option to avoid aborting when symmetry changes; 3 is supposed to be default according to documentation, but it seems that 0 (the default) is the only option that doesn't work from 0-4; also, it is interesting to note that all 4 options seem to work for test case with z-matrix input rather than xyz coords; cf. http://www.ccl.net/cgi-bin/ccl/message-new?2006+10+17+005 for original idea for solution
               "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3) IOP(4/21=2)", # use different SCF method; this addresses at least one case of failure for a C4H7J species
               "# pm3 opt=(verytight,calcfc,maxcyc=200) freq IOP(2/16=3) nosymm" , # try multiple different options (no gdiis, use calcfc, nosymm); 7/21/09: added maxcyc option to fix case of MPTBUKVAJYJXDE-UHFFFAOYAPmult3 (InChI=1/C4H10O5Si/c1-3-7-9-10(5,6)8-4-2/h4-5H,3H2,1-2H3/mult3) (file manually copied to speed things along)
               "# pm3 opt=(verytight,calcfc,maxcyc=200) freq=numerical IOP(2/16=3) nosymm", # numerical frequency keyword version of keyword #3; used to address GYFVJYRUZAKGFA-UHFFFAOYALmult3 (InChI=1/C6H14O6Si/c1-3-10-13(8,11-4-2)12-6-5-9-7/h6-7H,3-5H2,1-2H3/mult3) case; (none of the existing Gaussian or MOPAC combinations worked with it)
               "# pm3 opt=(verytight,gdiis,small) freq IOP(2/16=3)", #  somehow, this worked for problematic case of ZGAWAHRALACNPM-UHFFFAOYAF (InChI=1/C8H17O5Si/c1-3-11-14(10,12-4-2)13-8-5-7(9)6-8/h7-9H,3-6H2,1-2H3); (was otherwise giving l402 errors); even though I had a keyword that worked for this case, I manually copied the fixed log file to QMfiles folder to speed things along; note that there are a couple of very low frequencies (~5-6 cm^-1 for this case)
               "# pm3 opt=(verytight,nolinear,calcfc,small) freq IOP(2/16=3)", # used for troublesome C5H7J2 case (similar error to C5H7J below); calcfc is not necessary for this particular species, but it speeds convergence and probably makes it more robust for other species
               "# pm3 opt=(verytight,gdiis,maxcyc=200) freq=numerical IOP(2/16=3)", # use numerical frequencies; this takes a relatively long time, so should only be used as one of the last resorts; this seemed to address at least one case of failure for a C6H10JJ species; 7/15/09: maxcyc=200 added to address GVCMURUDAUQXEY-UHFFFAOYAVmult3 (InChI=1/C3H4O7Si/c1-2(9-6)10-11(7,8)3(4)5/h6-7H,1H2/mult3)...however, result was manually pasted in QMfiles folder to speed things along
               "# pm3 opt=tight freq IOP(2/16=3)", # this worked for problematic case of SZSSHFMXPBKYPR-UHFFFAOYAF (InChI=1/C7H15O5Si/c1-3-10-13(8,11-4-2)12-7-5-6-9-7/h7H,3-6H2,1-2H3) (otherwise, it had l402.exe errors); corrected log file was manually copied to QMfiles to speed things along; we could also add a freq=numerical version of this keyword combination for added robustness; UPDATE: see below
               "# pm3 opt=tight freq=numerical IOP(2/16=3)", # used for problematic case of CIKDVMUGTARZCK-UHFFFAOYAImult4 (InChI=1/C8H15O6Si/c1-4-12-15(10,13-5-2)14-7-6-11-8(7,3)9/h7H,3-6H2,1-2H3/mult4 (most other cases had l402.exe errors); corrected log file was manually copied to QMfiles to speed things along
               "# pm3 opt=(tight,nolinear,calcfc,small,maxcyc=200) freq IOP(2/16=3)", # similar to existing #5, but uses tight rather than verytight; used for ADMPQLGIEMRGAT-UHFFFAOYAUmult3 (InChI=1/C6H14O5Si/c1-4-9-12(8,10-5-2)11-6(3)7/h6-7H,3-5H2,1-2H3/mult3)
               "# pm3 opt freq IOP(2/16=3)", # use default (not verytight) convergence criteria; use this as last resort
               "# pm3 opt=(verytight,gdiis) freq=numerical IOP(2/16=3) IOP(4/21=200)", # to address problematic C10H14JJ case
               "# pm3 opt=(calcfc,verytight,newton,notrustupdate,small,maxcyc=100,maxstep=100) freq=(numerical,step=10) IOP(2/16=3) nosymm", # for very troublesome RRMZRNPRCUANER-UHFFFAOYAQ (InChI=1/C5H7/c1-3-5-4-2/h3H,1-2H3) case...there were troubles with negative frequencies, where I don't think they should have been; step size of numerical frequency was adjusted to give positive result; accuracy of result is questionable; it is possible that not all of these keywords are needed; note that for this and other nearly free rotor cases, I think heat capacity will be overestimated by R/2 (R vs. R/2) (but this is a separate issue)
               "# pm3 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm", #  for troublesome QDERTVAGQZYPHT-UHFFFAOYAHmult3(InChI=1/C6H14O4Si/c1-4-8-11(7,9-5-2)10-6-3/h4H,5-6H2,1-3H3/mult3); key aspects appear to be tight (rather than verytight) convergence criteria, no calculation of frequencies during optimization, use of numerical frequencies, and probably also the use of opt=small
               "# pm3 opt=(verytight,gdiis,calcall) IOP(2/16=3)", # used for troublesome C5H7J case; note that before fixing, I got errors like the following: "Incomplete coordinate system.  Try restarting with Geom=Check Guess=Read Opt=(ReadFC,NewRedundant) Incomplete coordinate system. Error termination via Lnk1e in l103.exe"; we could try to restart, but it is probably preferrable to have each keyword combination standalone; another keyword that may be helpful if additional problematic cases are encountered is opt=small; 6/9/09 note: originally, this had # pm3 opt=(verytight,gdiis,calcall) freq IOP(2/16=3)" (with freq keyword), but I discovered that in this case, there are two thermochemistry sections and cclib parses frequencies twice, giving twice the number of desired frequencies and hence produces incorrect thermo; this turned up on C5H6JJ isomer 
               "# pm3 opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm", # worked for troublesome ketene case: CCGKOQOJPYTBIH-UHFFFAOYAO (InChI=1/C2H2O/c1-2-3/h1H2) (could just increase number of iterations for similar keyword combination above (#6 at the time of this writing), allowing symmetry, but nosymm seemed to reduce # of iterations; I think one of nosymm or higher number of iterations would allow the similar keyword combination to converge; both are included here for robustness)
               "# pm3 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm", # added for case of ZWMVZWMBTVHPBS-UHFFFAOYAEmult3 (InChI=1/C4H4O2/c1-3-5-6-4-2/h1-2H2/mult3)
               "# pm3 opt=(calcall,small,maxcyc=100) IOP(2/16=3)", # used to address troublesome FILUFGAZMJGNEN-UHFFFAOYAImult3 case (InChI=1/C5H6/c1-3-5-4-2/h3H,1H2,2H3/mult3)
               ]

    @property
    def scriptAttempts(self):
        "The number of attempts with different script keywords"
        return len(self.keywords)

    @property
    def maxAttempts(self):
        "The total number of attempts to try"
        return 2 * len(self.keywords)

    def inputFileKeywords(self, attempt):
        """
        Return the top keywords for attempt number `attempt`.

        NB. `attempt`s begin at 1, not 0.
        """
        assert attempt <= self.maxAttempts
        if attempt > self.scriptAttempts:
            attempt -= self.scriptAttempts
        return self.keywords[attempt-1]

##########################################################################################

class GaussianTS(QMReaction, Gaussian):
    """
    A base Class for calculations of transition states using Gaussian. 

    Inherits from both :class:`QMReaction` and :class:`Gaussian`.
    """
    
    #: List of phrases to indicate success.
    #: ALL of these must be present in a successful job.
    successKeys = [
                   'Normal termination of Gaussian',
                   '******    1 imaginary frequencies (negative Signs) ******',
                  ]
    
    #: List of phrases that indicate failure
    #: NONE of these must be present in a succesful job.
    failureKeys = [
                   'ERROR TERMINATION',
                   'Error in internal coordinate system.',
                   ]
    
    def writeInputFile(self, attempt):
        """
        Using the :class:`Geometry` object, write the input file
        for the `attmept`th attempt.
        """
        chk_file = '%chk=' + os.path.join(self.settings.fileStore, self.uniqueID)
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "gjf")
        mol = openbabel.OBMol()
        
        obConversion.ReadFile(mol, self.geometry.getRefinedMolFilePath() )
        
        mol.SetTitle(self.uniqueID)
        obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
        input_string = obConversion.WriteString(mol)
        numProc = '%nprocshared=' + '4' # could be something that could be set in the qmSettings
        top_keys = self.keywords[0]
        title = ' ' + self.uniqueID
        with open(self.inputFilePath, 'w') as gaussianFile:
            gaussianFile.write(numProc)
            gaussianFile.write('\n')
            gaussianFile.write(chk_file)
            gaussianFile.write('\n')
            gaussianFile.write(top_keys)
            if attempt == 1:
                gaussianFile.write(input_string)
            else:
                gaussianFile.write('\n\n')
    
    def writeIRCFile(self):
        """
        Using the :class:`Geometry` object, write the input file for the 
        IRC calculation on the transition state. The geometry is taken 
        from the checkpoint file created during the geometry search.
        """
        chk_file = '%chk=' + os.path.join(self.settings.fileStore, self.uniqueID)
        numProc = '%nprocshared=' + '4' # could be something that could be set in the qmSettings
        top_keys = self.keywords[4]
        chrgMult = '1 2'
        with open(self.inputFilePath, 'w') as gaussianFile:
            gaussianFile.write(numProc)
            gaussianFile.write('\n')
            gaussianFile.write(chk_file)
            gaussianFile.write('\n')
            gaussianFile.write(top_keys)
            gaussianFile.write('\n\n')
            gaussianFile.write(chrgMult)
            gaussianFile.write('\n\n')
            
    def inputFileKeywords(self, attempt):
        """
        Return the top keywords.
        """
        raise NotImplementedError("Should be defined by subclass, eg. GaussianTSM062X")

    def generateQMKinetics(self):
        """
        Calculate the QM data and return a QMData object.
        """
        self.createGeometry()
        if self.verifyOutputFile():
            logging.info("Found a successful output file already; using that.")
        else:
            success = False
            for attempt in range(1, self.maxAttempts+1):
                self.writeInputFile(attempt)
                self.writeIRCFile()
                success = self.run()
                if success:
                    logging.info('Attempt {0} of {1} on species {2} succeeded.'.format(attempt, self.maxAttempts, self.molecule.toAugmentedInChI()))
                    break
            else:
                raise Exception('QM thermo calculation failed for {0}.'.format(self.molecule.toAugmentedInChI()))
        result = self.parse() # parsed in cclib
        return result
    
    def runIRC(self):
        self.testReady()
        # submits the input file to Gaussian
        process = Popen([self.executablePath, self.inputFilePath, self.outputFilePath])
        process.communicate()# necessary to wait for executable termination!
        
        return self.verifyIRCOutputFile()
    
    def verifyOutputFile(self):
        """
        Check's that an output file exists and was successful.
        
        Returns a boolean flag that states whether a successful GAUSSIAN simulation already exists for the molecule with the 
        given (augmented) InChI Key.
        
        The definition of finding a successful simulation is based on these criteria:
        1) finding an output file with the file name equal to the the reaction unique ID
        2) NOT finding any of the keywords that are denote a calculation failure
        3) finding all the keywords that denote a calculation success.
        
        If any of the above criteria is not matched, False will be returned and the procedures to start a new calculation 
        will be initiated.
        """
        if not os.path.exists(self.outputFilePath):
            logging.info("Output file {0} does not exist.".format(self.outputFilePath))
            return False
        
        # Initialize dictionary with "False"s 
        successKeysFound = dict([(key, False) for key in self.successKeys])
        failureKeysFound = dict([(key, False) for key in self.failureKeys])
        
        with open(self.outputFilePath) as outputFile:
            for line in outputFile:
                line = line.strip()
                
                for element in self.failureKeys: #search for failure keywords
                    if element in line:
                        logging.error("Gaussian output file contains the following error: {0}".format(element) )
                        failureKeysFound[element] = True
                    
                for element in self.successKeys: #search for success keywords
                    if element in line:
                        successKeysFound[element] = True
        
        if any(failureKeysFound.values()):
            if failureKeysFound['Error in internal coordinate system.']:
                return False, True
            else:
                return False, False
        
        # Check that ALL 'success' keywords were found in the file.
        if not all( successKeysFound.values() ):
            logging.error('Not all of the required keywords for success were found in the output file!')
            return False, False
        else:
            return True, False
    
    def verifyIRCOutputFile(self):
        """
        Check's that the resulting geometries of the path analysis match the reaction.
        """
        
        """
        Compares IRC geometries to input geometries.
        """
        if not os.path.exists(self.outputFilePath):
            logging.info("Output file {0} does not exist.".format(self.outputFilePath))
            return False
        
        # Initialize dictionary with "False"s 
        successKeysFound = dict([(key, False) for key in self.successKeys])
        
        pth1 = list()
        steps = list()
        with open(self.outputFilePath) as outputFile:
            for line in outputFile:
                line = line.strip()
                
                for element in self.failureKeys: #search for failure keywords
                    if element in line:
                        logging.error("Gaussian IRC output file contains the following error: {0}".format(element) )
                        return False
                    
                for element in self.successKeys: #search for success keywords
                    if element in line:
                        successKeysFound[element] = True
                
                if line.startswith('Point Number:'):
                    if int(line.split()[2]) > 0:
                        if int(line.split()[-1]) == 1:
                            ptNum = int(line.split()[2])
                            pth1.append(ptNum)
                        else:
                            pass
                elif line.startswith('# OF STEPS ='):
                    numStp = int(line.split()[-1])
                    steps.append(numStp)
        
        # Check that ALL 'success' keywords were found in the file.
        if not successKeysFound['Normal termination of Gaussian']:
            logging.error('Not all of the required keywords for success were found in the IRC output file!')
            return False
        # This indexes the coordinate to be used from the parsing
        elif steps == []:
            logging.error('No steps taken in the IRC calculation!')
            return False
        else:
            pth1End = sum(steps[:pth1[-1]])		
            # Compare the reactants and products
            ircParse = cclib.parser.Gaussian(self.outputFilePath)
            ircParse = ircParse.parse()
        
            atomnos = ircParse.atomnos
            atomcoords = ircParse.atomcoords
        
            # Convert the IRC geometries into RMG molecules
            # We don't know which is reactant or product, so take the two at the end of the
            # paths and compare to the reactants and products
            mol1 = Molecule()
            mol1.fromXYZ(atomnos, atomcoords[pth1End])
            mol2 = Molecule()
            mol2.fromXYZ(atomnos, atomcoords[-1])
            
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
                logging.warning("Didn't make expected reaction")
                return False
    
    def parseTS(self, labels):
    
        tsParse = cclib.parser.Gaussian(os.path.join(self.file_store_path, self.uniqueID + '.log'))
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
        
class GaussianTSM062X(GaussianTS):

    #: Keywords that will be added at the top of the qm input file
    keywords = [
               "# m062x/6-31+g(d,p) opt=(ts,calcall,tight,noeigentest)  int=ultrafine nosymm",
               "# m062x/6-31+g(d,p) opt=(ts,calcall,noeigentest) nosymm",
               "# m062x/6-31+g(d,p) irc=(calcall,report=read) geom=allcheck guess=check nosymm",
               ]
    """
    This needs some work, to determine options that are best used. Commented out the
    methods for now.
    """
    
    # @property
    # def scriptAttempts(self):
    #     "The number of attempts with different script keywords"
    #     return len(self.keywords)
    # 
    # @property
    # def maxAttempts(self):
    #     "The total number of attempts to try"
    #     return 2 * len(self.keywords)
    # 
    # def inputFileKeywords(self, attempt):
    #     """
    #     Return the top keywords for attempt number `attempt`.
    # 
    #     NB. `attempt`s begin at 1, not 0.
    #     """
    #     assert attempt <= self.maxAttempts
    #     if attempt > self.scriptAttempts:
    #         attempt -= self.scriptAttempts
    #     return self.keywords[attempt-1]

class GaussianTSB3LYP(GaussianTS):

    #: Keywords that will be added at the top of the qm input file
    keywords = [
               "# b3lyp/6-31+g(d,p) opt=(ts,calcall,tight,noeigentest) int=ultrafine nosymm",
               "# b3lyp/6-31+g(d,p) opt=(ts,calcall,tight,noeigentest,cartesian) int=ultrafine geom=allcheck guess=check nosymm",
               "# b3lyp/6-31+g(d,p) opt=(ts,calcall,noeigentest) nosymm",
               "# b3lyp/6-31+g(d,p) opt=(ts,calcall,noeigentest,cartesian) nosymm geom=allcheck guess=check nosymm",
               "# b3lyp/6-31+g(d,p) irc=(calcall,report=read) geom=allcheck guess=check nosymm",
               ]
    """
    This needs some work, to determine options that are best used. Commented out the
    methods for now.
    """

    # @property
    # def scriptAttempts(self):
    #     "The number of attempts with different script keywords"
    #     return len(self.keywords)
    # 
    # @property
    # def maxAttempts(self):
    #     "The total number of attempts to try"
    #     return 2 * len(self.keywords)
    # 
    # def inputFileKeywords(self, attempt):
    #     """
    #     Return the top keywords for attempt number `attempt`.
    # 
    #     NB. `attempt`s begin at 1, not 0.
    #     """
    #     assert attempt <= self.maxAttempts
    #     if attempt > self.scriptAttempts:
    #         attempt -= self.scriptAttempts
    #     return self.keywords[attempt-1]
