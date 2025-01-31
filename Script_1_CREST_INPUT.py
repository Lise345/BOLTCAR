#!/usr/bin/env python3
import sys
import re
import os
import shutil
import math
from unicodedata import digit
import numpy as np

#
# This script is intended to take as single argument a calculation in Gaussian 16:
#
#  - for TS a frequency calculation is the minimum criteria 
#
# From there it will ask the user the parameters required for CREST and DFT calculations
#
# It will then generate a submission file that will cascade through:
#
#  1 - Automated_RMSD (pruning + submitting Gaussian calculation) 
#  2 - Automated_IRC_Calculation (checking which TS has the correct frequency and running a forward and reverse IRC on them)
#  3 - Automated_IRC_Extractor (checking which IRC gave rise to a product and which gave rise to a complex and running geometry optimisations on them)
#  4 - TS_energies.py (optional, running the energy calculations on the TS for basis set extrapolation)
#  5 - Results_CREST.py (extracting results of the DFT calculations to an excel file and calculate Boltzmann population for each of the reactions) 
#

with open('./parameters.txt', 'r') as parameters:
    file_content = parameters.read()

    size_molecule = re.search(r'size_molecule\s*=\s*(\S+)', file_content)
    size_molecule = int(size_molecule.group(1))  # Adjust for zero-based indexing

# Get the file name from the command-line argument
file_name = sys.argv[1]
base_name = os.path.splitext(file_name)[0]

def lastgeometry(filename):
    with open(filename, "r") as readfile:
        lines = readfile.readlines()
        indices = []
        
        # Find all indices where "Standard orientation" appears
        for idx, line in enumerate(lines):
            if "Standard orientation" in line:
                indices.append(idx)
        
        # Get the last occurrence index
        if not indices:
            raise ValueError("No 'Standard orientation' found in the file.")
        last_index = indices[-1]
        
        # Extract geometry from lines after the last "Standard orientation"
        start = last_index + 5
        coord = []
        for i in range(start, start + size_molecule):
            strippedline = lines[i].split()
            number_list = [float(num) for num in strippedline]
            coord.append(number_list)
    
    return coord

def atomincoord(coord):
    for atom in coord:
        atom_number = atom[1]
        
        if atom_number in atomic_symbols:
            atom[1] = atomic_symbols[atom_number]
        else:
            print("Atom does not exist")
    return coord

# Dictionary mapping atomic numbers to element symbols 
atomic_symbols = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca",
    21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn",
    31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr",
    41: "Nb", 42: "Mo", 43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn",
    51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba", 57: "La", 58: "Ce", 59: "Pr", 60: "Nd",
    61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb",
    71: "Lu", 72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au", 80: "Hg",
    81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At", 86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th",
    91: "Pa", 92: "U", 93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk", 98: "Cf", 99: "Es", 100: "Fm",
    101: "Md", 102: "No", 103: "Lr", 104: "Rf", 105: "Db", 106: "Sg", 107: "Bh", 108: "Hs", 109: "Mt",
    110: "Ds", 111: "Rg", 112: "Cn", 113: "Nh", 114: "Fl", 115: "Mc", 116: "Lv", 117: "Ts", 118: "Og"
}

def convert_gjf_to_xyz(filename):
    coord=lastgeometry(filename)
    newcoord=atomincoord(coord)

    newfile="struc.xyz"
    if newfile in os.listdir():
        os.remove(newfile)
    xyzfile=open(newfile, "x")

    xyzfile.writelines(str(size_molecule) + '\n')
    xyzfile.writelines('\n')
    for atom in newcoord:
        xyzfile.writelines(atom[1]+' '+str(atom[3])+' '+str(atom[4])+' '+str(atom[5])+'\n')
    return 'Done'

# Check if calculation ended properly
with open(file_name, 'r') as f:
    file_content = f.read()

    if 'Normal termination of Gaussian 16' in file_content:
# Extract basic information about charge and multiplicity from the input file
        charge_match = re.search(r'Charge\s*=\s*(-?\d+)\s*', file_content)
        charge = int(charge_match.group(1))

        multiplicity_match = re.search(r'Multiplicity\s*=\s*(-?\d+)\s*', file_content)
        multiplicity = int(multiplicity_match.group(1))

        coordinates=lastgeometry(file_name)

# Create struc.xyz for CREST calculation
        with open('struc.xyz', 'w') as file:
            file.write(f'{size_molecule}\n')
            file.write('\n')
            file.write(f'{coordinates}')



# Extract imaginary frequency if it exists
        start_phrase = '******    1 imaginary frequencies (negative Signs) ******'
        end_phrase = 'Red. masses --'
        pattern = re.escape(start_phrase) + '(.*?)' + re.escape(end_phrase)
        matches = re.findall(pattern, file_content, re.DOTALL)

        if matches:
            frequencies = matches[-1]
            for line in frequencies.split('\n'):
                if 'Frequencies --' in line:
                    x = [float(num) for num in line.strip().split()[2:]]
                    print(x)
        else:
            x = None
            print(x)

# Stop the script if calculation did not converged
    else:
        print('Please provide a suitable input.')
        sys.exit()

convert_gjf_to_xyz(file_name)

# Get information from parameters.txt file
with open('parameters.txt', 'r') as parameters:
    file_content = parameters.read()

#######################################################################################################################################################################
#                                                                CREST parameters                                                                                     #
#######################################################################################################################################################################

# CREST version
    version_crest = re.search(r'CREST version(.+)', file_content) #mod FG 7/02/24: '=' removed

    if version_crest:
        version_crest = version_crest.group(1).strip() #mod FG 7/02/24: strip() added to remove trailing spaces
        print(version_crest)
        if version_crest == 'default': #mod FG 7/02/24
            version_crest = 'crest'
        elif version_crest == '3.0':
            version_crest = 'crest3'
        elif version_crest == 'continous release':
            version_crest = 'crest_continous_release'
        else:
            print('Please provide a suitable CREST version.')
            sys.exit()
    else:
        print('Unable to find CREST version information in the parameters file.')
        sys.exit()

# CREST solvent
    Solvent_CREST = re.search(r'CREST solvent(.+)', file_content)#mod FG 7/02/24: '=' removed

    if Solvent_CREST:
        Solvent_CREST = Solvent_CREST.group(1).strip()#mod FG 7/02/24: strip() added to remove trailing spaces

        if Solvent_CREST.lower() == 'none' :# mod FG 7/2/24: adding lowercase conversion and reduction of cases
            Solvent_CREST = ''
        elif re.search(r'--alpb', Solvent_CREST):
            pass
        else:
            print('Please provide a suitable CREST solvent.')
            sys.exit()
    else:
        print('Unable to find CREST solvent information in the parameters file.')
        sys.exit()

    NCI=re.search(r'NCI(.+)', file_content)

    if NCI:
        NCI=NCI.group(1).strip()
        if NCI.lower()=="none":
            NCI=""
        else:
            pass
    else:
        print('Unable to find NCI information in the parameters file.')
        sys.exit()
    

#######################################################################################################################################################################
#                                                                 RMSD parameters                                                                                     #
#######################################################################################################################################################################
        
# RMSD threshold
    RMSD = re.search(r'RMSD threshold(.+)' , file_content)

    if RMSD:
        RMSD = RMSD.group(1)

        try:
            RMSD_value = float(RMSD)
        except ValueError:
            print('Please provide a suitable RMSD threshold.')
            sys.exit()
    
    else:
        print('Unable to find RMSD threshold information in the parameters file.')
        sys.exit() 

#######################################################################################################################################################################
#                                                                  DFT parameters                                                                                     #
#######################################################################################################################################################################

# functionnal
    functional = re.search(r'Functional(.+)', file_content) #mod FG 7/02/24: '=' removed

    if functional:
        pass
    else:
        print('Unable to find functional information in the parameters file.')
        sys.exit()   

# base
    basis =  re.search(r'Basis(.+)', file_content)#mod FG 7/02/24: '=' removed
    if basis:
        pass
    else:
        print('Unable to find basis set information in the parameters file.')
        sys.exit()

# dispersion
    dispersion = re.search(r'Dispersion(.+)', file_content)#mod FG 7/02/24: '=' removed

    if dispersion:
        pass
    else:
        print('Unable to find dispersion information in the parameters file.')
        sys.exit() 

# solvent
    solvent = re.search(r'DFT solvent(.+)', file_content)#mod FG 7/02/24: '=' removed
    if solvent:
        pass
    else:
        print('Unable to find solvent information in the parameters file.')
        sys.exit()

# additional calculation
    AddCalc = re.search(r'Additional calculation (.+)', file_content)#mod FG 7/02/24: '=' removed
    if AddCalc:
        pass
    else:
        print('Unable to find additional calculation information in the parameters file.')
        sys.exit() 
        
# node for calculation
    node = re.search(r'Excluding nodes (.+)', file_content)#mod FG 7/02/24: '=' removed + string modified
    if node:
        pass
    else:
        print('Unable to find node excluded information in the parameters file.')
        sys.exit() 



# File with all parameters chosen
with open('parameters.txt', 'r') as parameters:
    file_content = parameters.read()

    pattern = re.compile(r'[-]+DFT parameters[-]+\n[-]+\n', re.MULTILINE)
    match = pattern.search(file_content)
    if x:
        modified_content = file_content[:match.end()] + '\n' + f'Experience name {base_name}' + '\n' + f'First frequency {int(x[0])}' + '\n' + f'Charge {charge}' + '\n' + f'Multiplicity {multiplicity}' + file_content[match.end():]
    else:
        modified_content = file_content[:match.end()] + '\n' + f'Experience name {base_name}' + '\n'+ f'First frequency none' + '\n' + f'Charge {charge}' + '\n' + f'Multiplicity {multiplicity}' + file_content[match.end():]
    
with open('parameters.txt', 'w') as file:
    file.write(modified_content)


# Construction of script.sub
multiplicityCREST = multiplicity - 1
with open('script_CREST.sub', 'w') as file:
    file.write('#!/bin/sh\n')
    file.write(f'#SBATCH --job-name=CREST_{base_name} \n')
    file.write('#SBATCH --nodes=1 \n')
    file.write('#SBATCH --ntasks=6 \n')
    file.write(f'#SBATCH --output=CREST_{base_name}.logfile \n')
    file.write('#SBATCH --time=48:00:00 \n')
    file.write('\n')
    file.write('cd ${SLURM_SUBMIT_DIR} \n')
    file.write('\n')
    file.write('#Loading modules \n')
    file.write('module load intel/2023a \n')
    file.write('module load AMS/2024.102-iimpi-2023a-intelmpi \n')
    file.write('module load sPyRMSD/0.8.0-foss-2023a\n')
    file.write('module load OpenBabel/3.1.1-gompi-2023a\n')
    file.write('module load SciPy-bundle/2023.07-gfbf-2023a\n')
    file.write('\n')
    file.write('#Launching calculation \n')
    file.write('\n')
    file.write(f'module load xtb/6.6.1-gfbf-2023a\n')
    file.write(f'module load CREST/2.12-gfbf-2023a\n')
    file.write(f'xtb struc.xyz --opt extreme --gfn 2 --input constraints.inp > xtb.out\n')
    file.write(f'crest xtbopt.xyz --T 6 --uhf {multiplicityCREST} --chrg {charge} {Solvent_CREST} {NCI} --cinp constraints.inp --subrmsd > CrestAnalysis.txt\n')	
    #file.write(f'crest struc.xyz --T 6 --uhf {multiplicityCREST} --chrg {charge} {Solvent_CREST} {NCI} --cinp constraints.inp --subrmsd > CrestAnalysis.txt\n')
    file.write(f'crest coord -cregen crest_conformers.xyz -ewin 30\n')
    #file.write(f'cp ../Script_1_RMSD_INPUT.py ./\n')
    file.write(f'./Script_1_RMSD_INPUT.py\n')
    file.write('\n')


# Construction of CREST folder
os.mkdir('CREST')
shutil.copy('struc.xyz', 'CREST/struc.xyz')
shutil.copy('parameters.txt', 'CREST/parameters.txt')
shutil.copy('script_CREST.sub', 'CREST/script_CREST.sub')
shutil.copy('Script_1_RMSD_INPUT.py', 'CREST/Script_1_RMSD_INPUT.py')
shutil.copy('constraints.inp', 'CREST/constraints.inp')

# Launch calculation
os.chdir('CREST')
os.system(f'sbatch script_CREST.sub')




