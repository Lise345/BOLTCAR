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
        ip.writelines("%nprocshared=4\n")
        ip.writelines("%mem=4GB\n")
        ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
        ip.writelines("# opt=(calcfc,ts,noeigentest) freq cc-pvdz empiricaldispersion=gd3 m062x\n")
        ip.writelines("\n")
        Title=filename[:-4]+" "+"cc-pVTZ_SP"+"\n"
        ip.writelines(Title)
        ip.writelines("\n")
        ip.writelines("0 1\n")
        
        for atom in lines[2:]:
            ip.writelines(atom)
        ip.writelines("\n")
        
        #Writing Link1 part for cc-pVTZ
        ip.writelines("--Link1--\n")
        ip.writelines("%nprocshared=4\n")
        ip.writelines("%mem=4GB\n")
        ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
        ip.writelines("# m062x cc-pvtz empiricaldispersion=gd3 Geom=Checkpoint \n")
        ip.writelines("\n")
        Title=filename[:-4]+" "+"cc-pVQT_SP"+"\n"
        ip.writelines(Title)
        ip.writelines("\n")
        ip.writelines("0 1\n")
        ip.writelines("\n")
        

	    #Writing Link1 part for cc-pVQZ
        ip.writelines("--Link1--\n")
        ip.writelines("%nprocshared=4\n")
        ip.writelines("%mem=4GB\n")
        ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
        ip.writelines("# m062x cc-pvqz empiricaldispersion=gd3 Geom=Checkpoint \n")
        ip.writelines("\n")
        Title=filename[:-4]+" "+"cc-pVQZ_SP"+"\n"
        ip.writelines(Title)
        ip.writelines("\n")
        ip.writelines("0 1\n")
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

    
