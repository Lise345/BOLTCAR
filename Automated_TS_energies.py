# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 14:41:14 2024

@author: killi
"""

import os
import sys
import subprocess
import shutil
import re
import numpy as np
import glob
import pandas as pd
import openpyxl

file_path=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "parameters.txt")
with open(file_path,'r') as parameters:
    file_content = parameters.read()
    
    functional = re.search(r'Functional (.+)', file_content)
    functional = functional.group(1).strip()

    basis_in= re.search(r'Basis (.+)', file_content)
    basis_in= basis.group(1).strip()
    if basis_in.lower()=='cbs':
        basis_2='cc-pvtz'
	basis_3='cc-pvqz'
    else:
        print("This calculation is only relevant for CBS extrapolation)
	sys.exit()

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

workbook = openpyxl.load_workbook("TS_analysis.xlsx")
sheet = workbook.active

data = []
for column in sheet.iter_cols(values_only=True):
    data.append(list(column))

print(data)

listofTS=[]
for caseofTS in data:
    if data.index(caseofTS)==2:
        for i in caseofTS:
            if caseofTS.index(i)>0 and i!=None:
                listofTS.append(i)

def SP_inputgenerator(xyzfile,filename):
    with open(xyzfile, 'r') as file:
        lines = file.readlines()
    with open(filename, 'x') as ip:
	#Writing part for cc-pVTZ
        ip.writelines("%nprocshared=4\n")
        ip.writelines("%mem=4GB\n")
        ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
        ip.writelines(f"# {functional} {basis_2} {dispersion} {solvent}\n")
        ip.writelines("\n")
        Title=filename[:-4]+" "+"cc-pVTZ_SP"+"\n"
        ip.writelines(Title)
        ip.writelines("\n")
        ip.writelines(f"{charge} {multiplicity}\n")
        
        for atom in lines[2:]:
            ip.writelines(atom)
        ip.writelines("\n")
        
        #Writing Link1 part for cc-pVQZ
        ip.writelines("--Link1--\n")
        ip.writelines("%nprocshared=4\n")
        ip.writelines("%mem=4GB\n")
        ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
        ip.writelines(f"# {functional} {basis_2} {dispersion} {solvent} Geom=Checkpoint \n")
        ip.writelines("\n")
        Title=filename[:-4]+" "+"cc-pVQZ_SP"+"\n"
        ip.writelines(Title)
        ip.writelines("\n")
        ip.writelines(f"{charge} {multiplicity}\n")
        ip.writelines("\n")
        ip.close()
    return print("Input generated for " + filename[:-4])

def launcher(xyzlist):
    for xyzfile in xyzlist:
        filename=xyzfile[:-4]+"_SP.gjf"
        with open(filename[:-4]+".sub","w") as gsub:
            gsub.write('#!/bin/sh\n')
            gsub.write(f'#SBATCH --job-name='+filename[:-4]+'\n')
            gsub.write('#SBATCH --ntasks=12\n')
            gsub.write(f'#SBATCH --output='+filename[:-4]+'.logfile\n')
            gsub.write('#SBATCH --time=01:00:00\n')
            gsub.write('\n')
            gsub.write('#Loading modules\n')
            gsub.write('module load intel/2023a\n')
            gsub.write('module load AMS/2024.102-iimpi-2023a-intelmpi\n')
            gsub.write('\n')
            gsub.write('export PATH=/vscmnt/brussel_pixiu_home/_user_brussel/105/vsc10536/bin/:$PATH\n')
            
            SP_inputgenerator(xyzfile,filename)
            
            gsub.write(f'sgx16 '+filename+' 47 \n')
            gsub.write('\n')
            gsub.close()
            
        os.system(f"sbatch "+filename[:-4]+".sub")
    return print("Launch of TS energies complete...")

launcher(listofTS)

    
