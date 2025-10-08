#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 14:16:23 2024

@author: lisevermeersch
"""


import os
import sys
import subprocess
import shutil
import re
import numpy as np
import glob
import pandas as pd

#-----------Loading parameters---------------

with open('./parameters.txt', 'r') as parameters:
    file_content = parameters.read()

    size_molecule = re.search(r'size_molecule\s*=\s*(\S+)', file_content)
    size_molecule = int(size_molecule.group(1)) - 1  # Adjust for zero-based indexing

    RMSD_threshold = re.search(r'RMSD threshold (.+)', file_content)
    RMSD_threshold = float(RMSD_threshold.group(1))

    Energy_threshold = re.search(r'Energy threshold (.+)', file_content)
    Energy_threshold = float(Energy_threshold.group(1))

    Energy_window = re.search(r'Energy window (.+)', file_content)
    Energy_window = float(Energy_window.group(1))

    B_threshold = re.search(r'B_threshold (.+)', file_content)
    B_threshold = float(B_threshold.group(1))

    rootdir = re.search(r'rootdir (.+)', file_content)
    rootdir = rootdir.group(1).strip().strip("'\"")

    binfolder = re.search(r'bin (.+)', file_content)
    binfolder = binfolder.group(1)

    #DFT parameters
    basis_in= re.search(r'Basis (.+)', file_content)
    basis_in= basis_in.group(1).strip()
    if basis_in.lower()=='cbs':
        basis_1='cc-pvdz'
    else:
        basis_1=basis_in

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

    functional = re.search(r'Functional (.+)', file_content)
    functional = functional.group(1).strip()

    time_irc = re.search(r'Time for IRC calcs\s+(\d+)', file_content)
    if time_irc:
        time_for_irc_calcs = int(time_irc.group(1))
        irc_time = f'{time_for_irc_calcs}:00:00'  # Format to HH:MM:SS
    else:
        irc_time = '25:00:00'  # Default value if not found
        print("Time for stationary calculations not found, defaulting to 25:00:00")

    # Extract the CC1_in value
    CC1_in_match = re.search(r'CC1_in\s*=\s*"(.*?)"', file_content)
    if CC1_in_match:
        CC1_in = list(map(int, CC1_in_match.group(1).split()))
    atom1=CC1_in[0]
    atom2=CC1_in[1]


#--------Check which TS's were found-------------

with open("CrestAnalysis.txt", "r") as crest:
    crestlines=crest.read()

    nconformers=re.search(r"number of unique conformers for further calc\s+(\d+)", crestlines)
    nconformers=nconformers.group(1)

    crestconformers=[]
    for i in range(1,int(nconformers)+1):
        name="crestconformer-"+str(i)
        crestconformers.append(name)

def compile_frequencies(lines):
    frequencies=[]
    for line in lines:
        if "Frequencies" in line:
            numbers = re.findall(r'-?\d+\.\d+', line)
            numbers = [float(num) for num in numbers]
            frequencies= frequencies+numbers
    return frequencies

def checkfrequency(filename):
    with open(filename, "r") as readfile:
        lines = readfile.readlines()
        last_line = lines[-1].rstrip()
        last_line_2 = lines[-2].rstrip()
    if "Normal termination" in last_line:
        frequencies=compile_frequencies(lines)
        nimag=0
        for frequency in frequencies:
            if frequency<0:
                nimag=nimag+1
        if nimag != 1:
            print("incorrect number of imag freq for "+filename)
            if nimag>1:
                incorrectTS.append(filename)
        elif frequencies[0]>-1000.0000 and frequencies[0]<-200.0000:
            correctTS.append(filename)
        else:
            incorrectTS.append(filename)
    else:
        errorterm.append(filename)

listfiles=[]
for file in os.listdir():
    if ".log" in file and ".logfile" not in file and 'IRC' not in file:
        listfiles.append(file)
            
incorrectTS= []
correctTS= []
errorterm=[]

for file in listfiles:
    checkfrequency(file)

def jobid(filename):
    match = re.findall(r"(\d+)", filename) #Finds all sequences of numbers
    if match:
        for group in match:
            if len(group)==4: #if sequence of numbers is 4 units long then return that as the jobid; unlikely that in name a sequence of 4 numbers would be added
                numbers=group
    return numbers

incorrectTS = sorted(incorrectTS, key=jobid)
correctTS = sorted(correctTS, key=jobid)
errorterm = sorted(errorterm, key=jobid)

print("The correct TS are")
print(correctTS)


print("The incorrect TS are")
print(incorrectTS)

print("The error TS are")
print(errorterm)

#------------Convert files to their xyz--------------
def lastgeometry(filename):
    with open(filename, "r") as readfile:
        lines = readfile.readlines()
        indices=[]
        for idx, line in enumerate(lines):  # Track index manually
            if "                         Standard orientation:                        " in line:
                indices.append(idx)
        last_index=indices[-1]
        print(last_index)
        
        start=last_index+5
        i=start
        coord=[]
        while (i-start)-1<size_molecule:
            strippedline=lines[i].split()
            number_list = [float(num) for num in strippedline]
            coord.append(strippedline)
            i=i+1
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

def atomincoord(coord):
    for atom in coord:
        atom_number = int(atom[1])
        
        if atom_number in atomic_symbols:
            atom[1] = atomic_symbols[atom_number]
        else:
            print("Atom does not exist")
    return coord


def convert_gjf_to_xyz(filename):
    coord=lastgeometry(filename)
    newcoord=atomincoord(coord)
    
    reduced_filename=filename[:].strip(".log")
    newfile=reduced_filename + ".xyz"
    if newfile in os.listdir():
        os.remove(newfile)
    xyzfile=open(newfile, "x")
    
    xyzfile.writelines(str(size_molecule) + '\n')
    xyzfile.writelines('\n')
    for atom in newcoord:
        xyzfile.writelines(atom[1]+' '+str(atom[3])+' '+str(atom[4])+' '+str(atom[5])+'\n')
    return 'Done'
        
for file in correctTS:
    convert_gjf_to_xyz(file)


#----------------Check Energies and Rotational constants for duplicate TS's-------------------

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

def xyzlistcleaner(filelist):
    print("RMSD cleaning in progress...")
    print("Total number of calculations="+str(len(filelist))) 
    toremovefromlist = []
    if len(filelist) > 1:
        for i, file1 in enumerate(filelist):
            if (i+1)<len(filelist):
                files_str=" ".join(filelist[(i+1):])
                
                result = subprocess.check_output(f'python -m spyrmsd -m {file1} {files_str}', shell=True) 
                
                RMSD_values = result.split()
                RMSD_values = [float(item.decode()) for item in RMSD_values]
                
                for idx, RMSD in enumerate(RMSD_values):
                    if 0<= RMSD <= float(RMSD_threshold):
                        file2=filelist[i+1+idx]
                        if file2 in toremovefromlist:
                            continue
                        else:
                            toremovefromlist.append(file2)
    print("Final removal list:", toremovefromlist)
    return toremovefromlist

def energyfinder(logfile):
    with open(logfile, 'r') as file:
        lines = file.readlines()
    indexlist=[]
    for line in lines:
        if "SCF Done" in line:
            indexlist.append(lines.index(line))
    lastindex=indexlist[len(indexlist)-1]
    linelist=lines[lastindex].split()
    for word in linelist:
        if word=="=":
            Energy=str(float(linelist[linelist.index(word)+1])*627.503)
    return Energy

def energycleaner(list):
    energies=[]
    for logfile in list:
        energy=energyfinder(logfile)
        energies.append(energy)

    minenergy=min(energies)
    maxenergy=max(energies)
    
    indices_to_remove = set()
    indices_converged=[]
    
    for i, energy_1 in enumerate(energies):
        for j, energy_2 in enumerate(energies):
            if j > i:
                energy_diff = abs(float(energy_1) - float(energy_2))
                if 0 <= energy_diff < Energy_threshold:
                    indices_converged.append([i,j])
                    indices_to_remove.add(j)        
                    print(f"Structures {list[i]} and {list[j]} have converged in terms of energy with a difference of {energy_diff} kcal/mol.")
                elif (float(minenergy)-float(energy_1))>Energy_window:
                    indices_to_remove.add(i)
                elif (float(minenergy)-float(energy_2))>Energy_window:
                    indices_to_remove.add(j)
                    
    cleanedlist = [logfile for k, logfile in enumerate(list) if k not in indices_to_remove]
    tooloworhighenergy = [logfile for k, logfile in enumerate(list) if k in indices_to_remove]
    converged_list = [[list[i], list[j]] for i, j in indices_converged]

    print("Structures that have converged in terms of energy: ")
    print(converged_list)

    print("Structures that have a distinguished energy: ")
    print(cleanedlist)
    
    print("Structures that have a too high or low energy: ")
    print(tooloworhighenergy)
    
    
    return [cleanedlist,converged_list,tooloworhighenergy]

def Bfinder(logfile):
    with open(logfile, "r") as file:
        lines=file.readlines()
    
    indices=[]
    for line in lines:
        if "Rotational constants (GHZ)" in line:
            index=lines.index(line)
            indices.append(index)
    
    Bline_index=indices[-1]
    Bline=lines[Bline_index]
    sBline=Bline.split()

    Bx=float(sBline[3])
    By=float(sBline[4])
    Bz=float(sBline[5])

    print("Rotational constants for "+logfile+" are "+str(Bx)+", "+str(By)+", "+str(Bz))
    
    return [Bx, By, Bz]

def Bcompiler(logfilelist):
    Bcomp=[]
    for logfile in logfilelist:
        Bees=Bfinder(logfile)
        Blist=[logfile,Bees]
        Bcomp.append(Blist)
    return Bcomp

def lowervalue(value1, value2):
    if value1< value2:
        return value1
    else:
        return value2

def Bcomparer(Bees1, Bees2):
    tot=0
    for i, B1 in enumerate(Bees1):
        for j, B2 in enumerate(Bees2):
            if i==j:
                lowest=lowervalue(float(B1),float(B2))
                diff=abs(float(B1)-float(B2))
                incr=(diff/lowest)*100
                sq_incr=incr**2
                tot=tot + sq_incr
    mean=tot/3
    RMSD=mean**(1/2)
    if RMSD<B_threshold:
        return True
    else:
        return False
                
def Bconverge(logfilelist):
    Bconv=[]
    B_ind = []
    for conformers_set in logfilelist:
        Bcomp=Bcompiler(conformers_set)
        log1 = Bcomp[0]
        i = 0
        for j, log2 in enumerate(Bcomp):
            if j>i and log1[0]!=log2[0]:
                #print(log1,log2)
                outcome=Bcomparer(log1[1], log2[1])
                #print(outcome)
                if outcome==True and log2[0] not in Bconv:
                    Bconv.append(log2[0])
                elif outcome == False and log2[0] not in B_ind and log2[0] not in Bconv:
                        B_ind.append(log2[0])
    B_ind = [el for el in B_ind if el not in Bconv]    
    return (B_ind, Bconv)

def cleaner(correctTS):
    #Energetic check
    print("Energy is being checked")
    energeticanalysis=energycleaner(correctTS)
    energycleaned=energeticanalysis[0]
    Econvergedlist=energeticanalysis[1]
    print("Energy has been evaluated")
    
    #Rotational Constant analysis
    print("Rotational constants are being checked")
    Banalysis=Bconverge(Econvergedlist)
    Bcleaned=Banalysis[0]
    Bconvergedlist=Banalysis[1]
    print("Rotational constants have been checked, BCleaned is")
    print(Bcleaned)
    print("Rotational constants have been checked, BConverged is")
    print(Bconvergedlist)
    
    IRClist=[]
    cleanedlogs = energycleaned + Bcleaned
    for log in cleanedlogs:
        xyzcorr=log[:-4]+".xyz"
        IRClist.append(xyzcorr)
    print("xyz files systems that are kept for IRC calculations are ")
    print(IRClist)
    return IRClist


IRClist=cleaner(correctTS)


#-------Save all information to an Excel file------------

def save_to_excel(incorrectTS, correctTS, IRClist, listfiles, crestconformers, errorterm, filename="TS_analysis.xlsx"):
    # Create a dictionary with lists to save
    data = {
        "CREST TS":pd.Series(crestconformers),
        "TS after optimization":pd.Series(listfiles),
        "Correct TS": pd.Series(correctTS),
        "IRC List": pd.Series(IRClist),
        "Incorrect TS": pd.Series(incorrectTS),
        "Error Termination": pd.Series(errorterm)
    }

    # Convert the dictionary into a DataFrame
    df = pd.DataFrame(dict([(k, pd.Series(v)) for k,v in data.items()]))

    # Save DataFrame to an Excel file
    df.to_excel(filename, index=False, engine='openpyxl')
    print(f"Data successfully saved to {filename}")

save_to_excel(incorrectTS, correctTS, IRClist, listfiles, crestconformers, errorterm)

print("Number of IRC calculations: "+str(len(IRClist)))

#------------IRC calculation input and submit script--------------------

def IRC_inputgenerator(xyzfile, filename, direction):
    with open(xyzfile, 'r') as file:
        lines = file.readlines()
    with open(filename, 'x') as ip:
        ip.writelines("%nprocshared=12\n")
        ip.writelines("%mem=12GB\n")
        ip.writelines(f"# irc=({direction},phase=({atom1},{atom2}),calcfc,maxpoints=100,recalc=3,Tight) {functional} {basis_1} {dispersion} {solvent}\n")
        ip.writelines("\n")
        Title=filename+" "+"IRC"+direction+"\n"
        ip.writelines(Title)
        ip.writelines("\n")
        ip.writelines(f"{charge} {multiplicity}\n")
        
        for atom in lines[2:]:
            ip.writelines(atom)
        ip.writelines("\n")
        ip.close()
    return print("IRC input generated for " + xyzfile[:-4])
        
def launcher(uplist,rootdir,binfolder):
    inp_file_job_ids = []
    
    n = 1
    
    for xyzfile in uplist:
        reduced_filename=xyzfile[:-4]+"_IRC"
        filename_forward=xyzfile[:-4]+"_IRCforward"+".gjf"
        filename_reverse=xyzfile[:-4]+"_IRCreverse"+".gjf"
        output_forward = filename_forward.replace(".gjf", ".log")
        output_reverse = filename_reverse.replace(".gjf", ".log")
        IRC_inputgenerator(xyzfile,filename_forward,"forward")
        IRC_inputgenerator(xyzfile,filename_reverse,"reverse")
        
        with open(reduced_filename+".sub","w") as gsub:
            gsub.write('#!/bin/sh\n')
            gsub.write(f'#SBATCH --job-name={reduced_filename}\n')
            gsub.write('#SBATCH --cpus-per-task=12\n')
            gsub.write(f'#SBATCH --output={reduced_filename}.logfile\n')
            gsub.write(f'#SBATCH --time={irc_time}\n')
            gsub.write('\n')
            gsub.write('# Loading modules\n')
            gsub.write('module load Gaussian/G16.A.03-intel-2022a\n')  # Adjust based on the available Gaussian module
            gsub.write('\n')
            gsub.write('# Setting up Gaussian environment\n')
            gsub.write('export GAUSS_SCRDIR=$VSC_SCRATCH_VO_USER/gauss_scrdir$SLURM_JOB_ID\n')  # Temporary directory for Gaussian scratch files
            gsub.write('mkdir -p $GAUSS_SCRDIR\n')
            gsub.write('#Launching calculation\n')
            gsub.write('export PATH={binfolder}:$PATH\n')
            gsub.write('dos2unix {filename_forward}\n')
            gsub.write('dos2unix {filename_reverse}\n')
            gsub.write(f'g16 < {filename_forward} > {output_forward} &\n')
            gsub.write(f'g16 < {filename_reverse} > {output_reverse} &\n')
            gsub.write('wait\n')
            gsub.write('rm -r ${VSC_SCRATCH_VO_USER:?}/gauss_scrdir${SLURM_JOB_ID:?}\n')
            gsub.write('\n')

    
        sbatch_command = f"sbatch {reduced_filename}.sub"
        
        
        result = subprocess.run(
            sbatch_command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        
        if result.returncode == 0:
        # Extract the job ID from the sbatch output
            job_id_match = re.search(r'(\d+)', result.stdout)
            if job_id_match:
                job_id = job_id_match.group(1)
                inp_file_job_ids.append(job_id)  # Collect job IDs for inp_file jobs
                n += 1  # Increment the counter
                print(f"Submitted job {n}")
            else:
                print(f"Failed to extract job ID for {sbatch_command}. Output: {result.stdout}")
        else:
            print(f"Failed to submit job for {sbatch_command}: {result.stderr}")
    
    if inp_file_job_ids:
        dependency_str = ":".join(inp_file_job_ids)
        extractor_script = os.path.join(rootdir,'3_StationaryPoints_calculator.sub')

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
    

launcher(IRClist,rootdir,binfolder)    

