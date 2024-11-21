#!/usr/bin/env python3
import os
import subprocess
import shutil
import re
import numpy as np
import glob

# This script is intended to be use after CREST calculation.
# The idea is to prune crest_conformers by using invariant RMSD with a threshold set to 0.5. https://github.com/RMeli/spyrmsd
# After pruning, calculation will be launch using DFT parameters define previously.

# Get information from parameters.txt file
with open('../parameters.txt', 'r') as parameters:
    file_content = parameters.read()

    # Extract the atoms as lists of integers
    molecule1_atoms = re.search(r'molecule1_atoms\s*=\s*(.+)', file_content)
    molecule1_atoms = list(map(int, molecule1_atoms.group(1).split()))

    molecule2_atoms = re.search(r'molecule2_atoms\s*=\s*(.+)', file_content)
    molecule2_atoms = list(map(int, molecule2_atoms.group(1).split()))

    # Extract atom indices and convert them to integers
    atom1_index = re.search(r'atom1\s*=\s*(\S+)', file_content)
    atom1_index = int(atom1_index.group(1)) - 1  # Adjust for zero-based indexing

    atom2_index = re.search(r'atom2\s*=\s*(\S+)', file_content)
    atom2_index = int(atom2_index.group(1)) - 1  # Adjust for zero-based indexing

    experience_number = re.search(r'Experience name (.+)', file_content)
    experience_number = experience_number.group(1)

    functional = re.search(r'Functional (.+)', file_content)
    functional = functional.group(1)

    base_1 =  re.search(r'First base (.+)', file_content)
    if base_1:
        base_1 = base_1.group(1)
        base_2 = re.search(r'Second base (.+)', file_content)
        base_2 = base_2.group(1)
    else:
        base_1 = re.search(r'Unique base (.+)', file_content)
        base_1 = base_1.group(1)
        base_2 = base_1

    pattern = re.compile('[\W_]+')
    base_name = pattern.sub('', base_1)

    dispersion = re.search(r'Dispersion (.+)', file_content)
    dispersion = dispersion.group(1)
    if dispersion == 'none' or dispersion == 'None':
        dispersion = ''

    solvent = re.search(r'DFT solvent (.+)', file_content)
    solvent = solvent.group(1)
    if solvent == 'none' or solvent == 'None':
        solvent = ''

    charge = re.search(r'Charge (-?\d+)', file_content)
    charge = charge.group(1)

    multiplicity = re.search(r'Multiplicity (-?\d+)', file_content)
    multiplicity = multiplicity.group(1)

    ref_freq = re.search(r'First frequency (.+)', file_content)
    if ref_freq.group(1).strip().lower() == 'none':
        ref_freq = ''
    else:
        ref_freq = int(ref_freq.group(1))
    node = re.search(r'Excluding nodes (.+)', file_content)
    node = node.group(1)
    if node == 'none' or node == 'None':
        node = ''
    
    RMSD_threshold = re.search(r'RMSD threshold (.+)', file_content)
    RMSD_threshold = RMSD_threshold.group(1)


def read_coordinates(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    start_index = 0
    for i, line in enumerate(lines):
        if '0 1' in line:
            start_index = i + 1
            break
    
    atoms = []
    for line in lines[start_index:]:
        parts = line.split()
        if len(parts) == 4:
            atoms.append([parts[0], float(parts[1]), float(parts[2]), float(parts[3])])
    return lines[:start_index], atoms

def write_coordinates(file_path, header, atoms):
    with open(file_path, 'w') as file:
        file.writelines(header)
        for atom in atoms:
            file.write(f"{atom[0]:<3} {atom[1]:>15.8f} {atom[2]:>15.8f} {atom[3]:>15.8f}\n")

def get_molecule_atoms(prompt):
    atoms = input(prompt).split()
    return [int(atom) - 1 for atom in atoms]

def get_atom(prompt):
    return int(input(prompt)) - 1

def process_file(file_path):
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    new_file_path = f"{base_name}_newgeom.xyz"
    
    header, atoms = read_coordinates(file_path)

    # Use extracted atom indices and coordinates
    atom1_coords = np.array(atoms[atom1_index][1:])
    atom2_coords = np.array(atoms[atom2_index][1:])
    vector = atom2_coords - atom1_coords
    distance = np.linalg.norm(vector)
    translation_vector = vector / distance * 0.4
    
    print(f"Translation vector: {translation_vector}")
    
    # Adjust atom indices for zero-based indexing
    for index in molecule1_atoms:
        index -= 1  # convert to zero-based indexing
        atoms[index][1:] = (np.array(atoms[index][1:]) + translation_vector).tolist()
    
    print("Translated Coordinates:")
    for index in molecule1_atoms:
        print(f"Atom {index+1}: {atoms[index]}")
    
    write_coordinates(new_file_path, header, atoms)
    print(f"The file has been updated successfully and saved as {new_file_path}.")

def XYZspliter():
    with open('crest_conformers.xyz', 'r') as rfile:
        lines = rfile.readlines()

    natoms = int(lines[0])
    ngeoms = len(lines) // (natoms + 2)

    for j in range(ngeoms):
        outname = f"xyzfilenum{j+1:04d}.xyz"
        with open(outname, "w") as ow:
            ow.write(str(natoms) + "\n \n")
            ow.writelines(lines[(j * (natoms + 2) + 2):((j + 1) * (natoms + 2))])

def xyzlistcleaner(filelist):
    toremovefromlist = []
    RMSD_value = []
    if len(filelist) > 1:

        result = subprocess.check_output(f'python -m spyrmsd -m {filelist[0]} xyzfilenum*', shell=True)

        RMSD_value = result.split()
        RMSD_value = [item.decode() for item in RMSD_value]
        RMSD_value.extend(RMSD_value)

        for i in range(1,len(filelist)):
            if 0 < float(RMSD_value[i]) <= float(RMSD_threshold):
                toremovefromlist.append(filelist[i])

    return toremovefromlist

with open('CrestAnalysis.txt', 'r') as readfile:
    lines = readfile.readlines()
    last_line = lines[-1].rstrip()
    last_line_2 = lines[-2].rstrip()

if ' CREST terminated normally.' in last_line:
    if ref_freq:
        os.mkdir('./TS-CREST')
        shutil.copy('crest_conformers.xyz', './TS-CREST/crest_conformers.xyz') 
        with open('../log.txt', 'a') as log:
            log.write('CREST completed\n')

        os.chdir('./TS-CREST')
    else:
        os.mkdir('./Geo-CREST')
        shutil.copy('crest_conformers.xyz', './Geo-CREST/crest_conformers.xyz') 
        with open('../log.txt', 'a') as log:
            log.write('CREST completed\n')

        os.chdir('./Geo-CREST')

    for file_path in glob.glob('./*.xyz'):
        print(f"Processing file: {file_path}")
        process_file(file_path)

    XYZspliter()

    dirs = os.listdir()
    xyzlist = [ele for ele in dirs if ".xyz" in ele and "num" in ele]
    ordxyzlist = sorted(xyzlist)

    n = 1
    listn = ordxyzlist[:]
    print(listn)

    while len(listn) > 0:
        match = re.search(r'\d+', listn[0]).group()
        fout = open(f"{base_name}_{experience_number}-{match}.inp", "w+")
        with open(listn[0], 'r') as fp:
            text = fp.read().splitlines(True)[2:]
            if ref_freq:
                fout.writelines("%nprocshared=12\n")
                fout.writelines("%mem=5GB\n")
                fout.writelines(f"# opt=(calcfc,ts,noeigen) freq {functional} {base_1} {dispersion}\n")
                fout.writelines("\n")
                fout.writelines(f"H2\n")
                fout.writelines(f"\n")
                fout.writelines(f"{charge} {multiplicity}\n")
            else:
                fout.writelines("%nprocshared=12\n")
                fout.writelines("%mem=5GB\n")
                fout.writelines(f"# opt {functional} {base_1} {dispersion} \n")
                fout.writelines("\n")
                fout.writelines(f"H2\n")
                fout.writelines(f"\n")
                fout.writelines(f"{charge} {multiplicity}\n")
            for item in text:
                fout.writelines(''.join(map(str, item)))
            fout.writelines('\n')
        fp.close()

        removals = xyzlistcleaner(listn)
        for file_to_remove in removals:
            os.remove(file_to_remove)
        listn = [item for item in listn if item not in removals]
        os.remove(listn[0])
        listn = listn[1:]
        
        with open(f'{base_name}_{experience_number}-{match}g16.sub', 'w') as gsub:
            gsub.write('#!/bin/sh\n')
            gsub.write(f'#SBATCH --job-name={base_name}_{experience_number}-{match}\n')
            gsub.write('#SBATCH --ntasks=12\n')
            gsub.write(f'#SBATCH --output={base_name}_{experience_number}-{match}.logfile\n')
            gsub.write('#SBATCH --time=05:00:00\n')
            gsub.write('\n')
            gsub.write('#Loading modules\n')
            gsub.write('module load intel/2023a\n')
            gsub.write('module load AMS/2024.102-iimpi-2023a-intelmpi\n')
            gsub.write('\n')
            gsub.write('#Launching calculation\n')
            gsub.write('export PATH=/vscmnt/brussel_pixiu_home/_user_brussel/105/vsc10536/bin/:$PATH\n')
            gsub.write(f'sgx16 {base_name}_{experience_number}-{match}.inp 5 \n')
            gsub.write('\n')

        os.system(f"sbatch {base_name}_{experience_number}-{match}g16.sub")

        n += 1
