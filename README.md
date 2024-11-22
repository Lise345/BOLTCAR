# FASTCAR_Adapted
This code serves to analyse reactions of structures with flexible side chains. By using CREST to study the transition state structure (TS), all different kind of conformers can be studied. Our adaptation of the FASTCAR method automatically performs the CREST method and subsequently optimizes the corresponding transition state structures. To make sure these don't correspond to the same structure, the energy of the TS is analyzed. If two TS have a similar energy, the rotational constants will be evaluated to check whether the found TS are indeed the same structures. If so, these TS's are discarted. Afterwards the IRC of the structure is calculated. Then the complex and product structures are optimized and depending on the settings in the parameters.txt file the energies are extrapolated to an infinite basis set or not. Finally, the results are summarized in the adapted_FASTCAR_results.xlsx file, where the Boltzmann distribution is also used to determine the prevalence of each reaction. 
# Startgeom.out
The startgeom.out contains the output file of a Gaussian 16 calculation, where the initial TS was optimized. This will serve as the starting point of your calculation.
# Parameters.txt
This file will contain parameters used throughout the code:
CREST version --> indicate the CREST version
CREST solvent --> indicate the CREST solvent
RMSD threshold --> classically equal to 0.5
Functional --> used functional for Gaussian optimization
Unique base --> basis set for Gaussian optimization
DFT solvent --> indicate solvent if necessary
Additional calculation --> ?
Excluding nodes --> 
molecule1_atoms = --> list the atoms in molecule 1, space separated
molecule2_atoms = --> list the atoms in molecule 2, space separated
atom1 = --> indicate one atom that participates in the bond formation, coming from molecule 1
atom2 = --> indicate the atom that will be bonded to atom 1
size_molecule = --> indicate the number of atoms in your molecule
CC1_in = --> indicate atom 1 and atom 2, e.g. "1 2"
# Constraints.inp



