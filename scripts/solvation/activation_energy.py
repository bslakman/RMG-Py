# Script written by Pierre, and modified by Belinda.

from cclib import parser
import argparse
import collections
import os

solventList = ['', '_n-octane']#["", "_water", "_n-octane", "_benzene", "_pyridine", "_tetrahydrofuran", "_dichloromethane", "_acetonitrile", "_dimethylsulfoxide"]

clParser = argparse.ArgumentParser(description="""
Given a reaction family, will get TS energies from Gaussian output files in
different solvents
""")
clParser.add_argument("-f", "--family", default="H_Abstraction", help="Name of the family")
#clParser.add_argument("-r", "--reactants", nargs='+', help="list of reactants (strings)")
#clParser.add_argument("-o", "--optional", default = "", nargs='+', help="""\
#optional stuff that tells you about the reacting site like CH or prim...
#""")
args = clParser.parse_args()
directory = os.path.join("/Users/belinda/Gaussian/SMD/", args.family)

# Go through all of the reactions for that family
for rxn_folder in os.listdir(directory):
    rxn_directory = os.path.join(directory, rxn_folder)

    if '+' in rxn_folder: # bimolecular
        r1 = rxn_folder.split('_')[0].split('+')[0]
        r2 = rxn_folder.split('_')[0].split('+')[1]
    else:
        r1 = rxn_folder.split('_')[0] # unimolecular
        r2 = None

    # Do for all the solvents in the list.
    for solvent in solventList:
        reactantOutput = os.path.join(rxn_directory, r1 + solvent + ".log")
        if r2:
            reactant2Output = os.path.join(rxn_directory, r2 + solvent + ".log")
        else: reactant2Output = None
        tsOutput = os.path.join(rxn_directory, "ts" + solvent + ".log")

        outputDataFile = os.path.join(rxn_directory, "output" + solvent+ ".txt")

        if os.path.exists(reactantOutput) and os.path.exists(tsOutput) and not os.path.exists(outputDataFile):
            rParse = parser.Gaussian(reactantOutput)
            tsParse = parser.Gaussian(tsOutput)

            try:
                rParse = rParse.parse()
            except AssertionError:
                continue
            try:
                tsParse = tsParse.parse()
            except AssertionError:
                continue

            # In Hartrees
            reactantE = rParse.scfenergies[-1]/27.2113845
            tsE = tsParse.scfenergies[-1]/27.2113845
            tsVib = tsParse.vibfreqs[0]

            if reactant2Output is not None:
                r2Parse = parser.Gaussian(reactant2Output)
                try:
                    r2Parse = r2Parse.parse()
                except AssertionError:
                    continue
                reactant2E = r2Parse.scfenergies[-1]/27.2113845
            else:
                reactant2E = 0.0

            Ea = (tsE - reactantE - reactant2E) * 2600
            if solvent is "":
                gasEa = Ea
            try:
                diffEa = Ea - gasEa
            except NameError:
                continue

            rString = 'Reactant energy = ' + str(reactantE)
            r2String = 'Reactant 2 energy = ' + str(reactant2E)
            tEnergy = 'TS energy       = ' + str(tsE)
            EaString = 'Activation energy (in kJ/mol)     = ' + str(Ea)
            tVib    = 'TS vib          = ' + str(tsVib)
            diffString = 'Difference in activation energy from gas phase (in kJ/mol)   = ' + str(diffEa)

            with open(outputDataFile, 'w') as parseFile:
                parseFile.write('The energies of the species in Hartree are:')
                parseFile.write('\n')
                parseFile.write(rString)
                parseFile.write('\n')
                parseFile.write(r2String)
                parseFile.write('\n')
                parseFile.write(tEnergy)
                parseFile.write('\n')
                parseFile.write(tVib)
                parseFile.write('\n')
                parseFile.write(EaString)
                parseFile.write('\n')
                parseFile.write(diffString)
