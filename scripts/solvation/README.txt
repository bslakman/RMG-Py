1) solvation_energies.py: get TS energies in solvent (on discovery)
2) reactant_energies.py: get reactant energies in solvent; also copies over gas phase files from Pierre’s scratch (on discovery)
3) Copy data folder onto local computer (not necessary?. we have cclib on discovery)
4) activation_energy.py: get output_solvent.txt files for every reaction (local)
5) SMD.py: get all the reactions and Ea data and dictionaries into RMG-database folders (local)
6) Train the groups
     -copy over groups.py to solvationGroups.py for now;
     -remove template, recipe and forbidden groups;
     -replace “kinetics” with “correction"
     -run solvationGroups.py (in RMG-database folder)
