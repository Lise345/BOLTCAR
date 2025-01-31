# FASTCAR_Adapted

This code serves to analyse reactions of structures with flexible side chains. By using CREST to study the transition state structure (TS), all different kind of conformers can be studied. Our adaptation of the FASTCAR method automatically performs the CREST method and subsequently optimizes the corresponding transition state structures. To make sure these don't correspond to the same structure, the energy of the TS is analyzed. If two TS have a similar energy, the rotational constants will be evaluated to check whether the found TS are indeed the same structures. If so, these TS's are discarted. Afterwards the IRC of the structure is calculated. Then the complex and product structures are optimized and depending on the settings in the parameters.txt file the energies are extrapolated to an infinite basis set or not. Finally, the results are summarized in the adapted_FASTCAR_results.xlsx file, where the Boltzmann distribution is also used to determine the prevalence of each reaction. 

# IMPORTANT

Check all submit scripts (.sub) as you will have to make changes to select your path, as well as the parameters.txt file.

# NAMEOFTHEFILE.out

The startgeom.out contains the output file of a Gaussian 16 calculation, where the initial TS was optimized. This will serve as the starting point of your calculation.

# Constraints.inp

This file is a standard constraints file for a CREST calculation. In our example we fixed two bonds for which the distances are set to those in the Diels-Alder transition state. Make sure the constraints are not too loose, or you will find unrealistic results.

# Parameters.txt

The Parameters file is used throughout the code and will set several variables that will be used.

--------------------------------------------------
-----------------CREST parameters-----------------
--------------------------------------------------

CREST version --> select CREST version  

CREST solvent --> type solvent if necessary  

--------------------------------------------------
------------------RMSD parameters-----------------
--------------------------------------------------

RMSD threshold  0.5 --> maximum RMSD threshold used to see if some TS's have converged
Energy threshold 0.05 --> maximum energy threshold used to see if some TS's have converged (in kcal/mol)
Energy window 10 --> maximum energy window used to see if some TS's have converged (in kcal/mol)
B_threshold 2 --> How to define this?

--------------------------------------------------
------------------DFT parameters------------------
--------------------------------------------------

Functional  m062x --> Type functional of choice
Dispersion  empiricaldispersion=gd3 --> Type dispersion of choice or None if none should be used
Basis cbs --> Type basis set of choice, "cbs" will do an infinite basis set extrapolation using cc-pvdz/tz/qz basis sets
DFT solvent  none --> Type solvent of choice or none
Energy of separate reagents -1200.27851165566 --> Type energy of separate reagents (will only be used if separate reagents are not calculated (TODO))
Charge 0 --> Type Charge of system
Multiplicity 1 --> Type Multiplicity of system
Time for IRC calcs 40 --> Type time required for IRC calcs
Time for stationary calcs 40 --> Type time required for Stationary calcs
Time for TS calcs 40 --> Type time required for TS calcs

--------------------------------------------------
---------------Geometry information---------------
--------------------------------------------------

molecule1_atoms = 3 8 9 12 13 14 17 18 19 20 21 22 23 24 42 43 44 45 --> Type atoms in molecule 1 (for calc of energy of separate reagents)
molecule2_atoms = 1 2 4 5 6 7 10 11 15 16 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 --> Type atoms in molecule 2 (for calc of energy of separate reagents)
atom1 = 1 --> Type atom that is part of C-C bond in molecule 1
atom2 = 12 --> Type atom that is part of C-C bond in molecule 2, connected to atom 1
size_molecule = 45 --> Type number of atoms in the molecule
CC1_in = "1 12" --> Type atoms that are part of C-C bond

--------------------------------------------------
---------------------Used path--------------------
--------------------------------------------------

rootdir /rhea/scratch/brussel/105/vsc10536/lise/13_Ondemand/Fulvenes/Fullstructure/P_CyFv/BOLTCAR --> Type folder in which you will perform BOLTCAR
bin /vscmnt/brussel_pixiu_home/_user_brussel/105/vsc10536/bin/ --> Type path to your bin folder



