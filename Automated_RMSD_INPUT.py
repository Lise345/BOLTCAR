#!/usr/bin/env python3
import os
import subprocess
import shutil
import re
import numpy as np
import glob
import sys

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
    functional = functional.group(1).strip()

    basis_1 =  re.search(r'Basis (.+)', file_content)
    basis_1 = basis_1.group(1).strip()
    if basis_1.lower()=="cbs":
        basis_1="cc-pvdz"

    pattern = re.compile('[\W_]+')
    base_name = pattern.sub('', basis_1)

    dispersion = re.search(r'Dispersion (.+)', file_content)
    dispersion = dispersion.group(1).strip()
    if dispersion == 'none' or dispersion == 'None':
        dispersion = ''

    solvent = re.search(r'DFT solvent (.+)', file_content)
    solvent = solvent.group(1).strip()
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
        if ref_freq >= 0:
            print("The structure that was provided is not a transition state. Please provide one.")
            sys.exit()
    
    node = re.search(r'Excluding nodes (.+)', file_content)
    node = node.group(1)
    if node == 'none' or node == 'None':
        node = ''
    
    RMSD_threshold = re.search(r'RMSD threshold (.+)', file_content)
    RMSD_threshold = RMSD_threshold.group(1)

    rootdir = re.search(r'rootdir (.+)', file_content)
    rootdir = rootdir.group(1).strip().strip("'\"")

    binfolder = re.search(r'bin (.+)', file_content)
    binfolder = binfolder.group(1)


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
    os.mkdir('./TS-CREST')
    shutil.copy('crest_conformers.xyz', './TS-CREST/crest_conformers.xyz') 
    with open('../log.txt', 'a') as log:
        log.write('CREST completed\n')

    os.chdir('./TS-CREST')

    for file_path in glob.glob('./*.xyz'):
        print(f"Processing file: {file_path}")
        process_file(file_path)

    XYZspliter()

    dirs = os.listdir()
    xyzlist = [ele for ele in dirs if ".xyz" in ele and "num" in ele]
    ordxyzlist = sorted(xyzlist)

    n = 1
    MAX_JOBS = 10
    MAX_SUB_FILES = 50  # Maximum number of .sub files allowed
    listn = ordxyzlist[:]
    print(listn)
    
    inp_file_job_ids = []
    
    while listn and n <= MAX_SUB_FILES:  # Stop if the limit is reached
        match = re.search(r'\d+', listn[0]).group()
        fout = open(f"input_{experience_number}-{match}.inp", "w+")
        with open(listn[0], 'r') as fp:
            text = fp.read().splitlines(True)[2:]
            fout.writelines("%nprocshared=12\n")
            fout.writelines("%mem=16GB\n")
            fout.writelines(f"# opt=(calcfc,ts,noeigen) freq {functional} {basis_1} {dispersion} {solvent}\n")
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
            gsub.write('# Loading modules\n')
            gsub.write('module load Gaussian/G16.A.03-intel-2022a\n')
            gsub.write('\n')
            gsub.write('# Setting up Gaussian environment\n')
            gsub.write('export GAUSS_SCRDIR=$TMPDIR\n')  # Temporary directory for Gaussian scratch files
            gsub.write('mkdir -p $GAUSS_SCRDIR\n')
            gsub.write('#Launching calculation\n')
            gsub.write('export PATH={binfolder}:$PATH\n')
            gsub.write('dos2unix input_{experience_number}-{match}.inp\n')
            gsub.write(f'g16 < input_{experience_number}-{match}.inp > {base_name}_{experience_number}-{match}.log\n')
            gsub.write('\n')
    
        sbatch_command = f"sbatch {base_name}_{experience_number}-{match}g16.sub"
        
        result = subprocess.run(
            sbatch_command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        
        if result.returncode == 0:
            job_id_match = re.search(r'(\d+)', result.stdout)
            if job_id_match:
                job_id = job_id_match.group(1)
                inp_file_job_ids.append(job_id)
                n += 1
                print(f"Submitted job {n}/{MAX_JOBS} for file")
            else:
                print(f"Failed to extract job ID for {sbatch_command}. Output: {result.stdout}")
        else:
            print(f"Failed to submit job for {sbatch_command}: {result.stderr}")
    
    if inp_file_job_ids:
        dependency_str = ":".join(inp_file_job_ids)
        extractor_script = os.path.join(rootdir, '2_IRC_calculator.sub')
    
        if not os.path.exists(extractor_script):
            raise FileNotFoundError(f"Extractor script not found: {extractor_script}")
    
        dependency_command = [
            "sbatch",
            f"--dependency=afterany:{dependency_str}",
            extractor_script
        ]
    
        result = subprocess.run(dependency_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if result.returncode == 0:
            print(f"Extractor job submitted successfully: {result.stdout}")
        else:
            print(f"Failed to submit extractor job: {result.stderr}")
    else:
        print("No jobs were submitted, skipping dependency job submission.")

        

       
