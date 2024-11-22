import os
import math
import sys
import re

with open('../parameters.txt', 'r') as parameters:
    file_content = parameters.read()
    
    functional = re.search(r'Functional (.+)', file_content)
    functional = functional.group(1)

    basis_in= re.search(r'Basis (.+)', file_content)
    basis_in= basis.group(1)
    if basis_in.lower()=='cbs':
        basis_1='cc-pvdz'
        basis_2='cc-pvtz'
        basis_3='cc-pvqz'
    else:
        basis_1=basis_in

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
    
    # Extract the size_molecule value
    size_molecule_match = re.search(r'size_molecule\s*=\s*(\d+)', file_content)
    if size_molecule_match:
        size_molecule = int(size_molecule_match.group(1))

    # Extract the CC1_in value
    CC1_in_match = re.search(r'CC1_in\s*=\s*"(.*?)"', file_content)
    if CC1_in_match:
        CC1_in = list(map(int, CC1_in_match.group(1).split()))

CC1_out=CC1_in
CC1=CC1_in



def geometryextractor(logfile):
    with open(logfile, 'r') as file:
        lines = file.readlines()

    convergence_indices = []
    for i, line in enumerate(lines):
        if "Input orientation" in line:
            convergence_indices.append(i)

    conv_geom=convergence_indices[-1]+5

    print("geometry found ...")
    
    geometry= []
    start=conv_geom
    index=conv_geom
    while (index-start) < size_molecule:
        geometry.append(lines[index])
        index=index+1
    return geometry

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

def geometryconverter(geometry):
    updated_geometry=[]
    for atom in geometry:
        updated_atom=atom.split()
        atom_number = int(updated_atom[1])
        if atom_number in atomic_symbols:
            updated_atom[1] = atomic_symbols[atom_number]
        else:
            print("Atom does not exist")
        del updated_atom[2]
        updated_geometry.append(updated_atom)
    print(updated_geometry)
    return updated_geometry

def inputgenerator(geometry, filename):
    with open(filename, 'x') as ip:
        ip.writelines("%nprocshared=8\n")
        ip.writelines("%mem=16GB\n")
        ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
        ip.writelines(f"# opt=calcfc freq {functional} {basis_1} {dispersion} {solvent}\n")
        ip.writelines("\n")
        Title=filename[:-4]+" "+"optfreq"+"\n"
        ip.writelines(Title)
        ip.writelines("\n")
        ip.writelines(f"{charge} {multiplicity}\n")
        
        for atom in geometry:
            atom_line=atom[1]+" "+atom[2]+" "+atom[3]+" "+atom[4]+"\n"
            ip.writelines(atom_line)
        ip.writelines("\n")

        if basis_in.lower()=="cbs":
        #Writing Link1 part for cc-pVTZ
            ip.writelines("--Link1--\n")
            ip.writelines("%nprocshared=4\n")
            ip.writelines("%mem=4GB\n")
            ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
            ip.writelines(f"# {functional} {basis_2} {dispersion} {solvent} Geom=Checkpoint\n")
            ip.writelines("\n")
            Title=filename[:-4]+" "+"E_ccpvtz"+"\n"
            ip.writelines(Title)
            ip.writelines("\n")
            ip.writelines(f"{charge} {multiplicity}\n")
            ip.writelines("\n")
            
            #Writing Link1 part for cc-pVQZ
            ip.writelines("--Link1--\n")
            ip.writelines("%nprocshared=8\n")
            ip.writelines("%mem=4GB\n")
            ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
            ip.writelines(f"# {functional} {basis_3} {dispersion} {solvent} Geom=Checkpoint\n")
            Title=filename[:-4]+" "+"E_ccpvqz"+"\n"
            ip.writelines("\n")
            ip.writelines(Title)
            ip.writelines("\n")
            ip.writelines(f"{charge} {multiplicity}\n")
            ip.writelines("\n")
            ip.close()
    return print("Input generated for " + filename[:-4])

def atomwithfloats(atom):
    updated_atom=[]
    for el in atom:
        if atom.index(el)==0 or atom.index(el)==1:
            updated_atom.append(el)
        else:
            up_el=float(el)
            updated_atom.append(up_el)
    return updated_atom

def distance(geometry, coordinates):
    atom1=atomwithfloats(geometry[coordinates[0]])
    atom2=atomwithfloats(geometry[coordinates[1]])
    r_sq=abs(((atom1[2]-atom2[2])**2)+((atom1[3]-atom2[3])**2)+((atom1[4]-atom2[4])**2))
    
    r=math.sqrt(r_sq)
    return r


def launcher(logfilelist):
    for i, logfile1 in enumerate(logfilelist):
        number=logfile1[17:21]
        for j, logfile2 in enumerate(logfilelist):
            if number in logfile2 and logfile1 != logfile2 and j>i:
                geometry1=geometryextractor(logfile1)
                updated_geometry1=geometryconverter(geometry1)
                distance1_CC1=distance(updated_geometry1, CC1)
                
                geometry2=geometryextractor(logfile2)
                updated_geometry2=geometryconverter(geometry2)
                distance2_CC1=distance(updated_geometry2, CC1)
                
                reduced_filename=logfile1[:-4]+"_optE"
                with open(reduced_filename+".sub","w") as gsub:
                    gsub.write('#!/bin/sh\n')
                    gsub.write(f'#SBATCH --job-name='+reduced_filename+'\n')
                    gsub.write('#SBATCH --ntasks=12\n')
                    gsub.write(f'#SBATCH --output='+reduced_filename+'.logfile\n')
                    gsub.write('#SBATCH --time=01:00:00\n')
                    gsub.write('\n')
                    gsub.write('#Loading modules\n')
                    gsub.write('module load intel/2023a\n')
                    gsub.write('module load AMS/2024.102-iimpi-2023a-intelmpi\n')
                    gsub.write('\n')
                    gsub.write('export PATH=/vscmnt/brussel_pixiu_home/_user_brussel/105/vsc10536/bin/:$PATH\n')
                
                    if distance1_CC1 > distance2_CC1:
                        filename1=logfile1[:-4]+"_Complex.gjf"
                        filename2=logfile2[:-4]+"_Product.gjf"
                    
                        inputgenerator(updated_geometry1, filename1)
                        inputgenerator(updated_geometry2, filename2)
                        
                        gsub.write(f'sgx16 '+filename1+' 119 \n')
                        gsub.write(f'sgx16 '+filename2+' 119 \n')
                        gsub.write('\n')
                        
                    else:
                        filename2=logfile2[:-4]+"_Complex.gjf"
                        filename1=logfile1[:-4]+"_Product.gjf"
                    
                        inputgenerator(updated_geometry1, filename1)
                        inputgenerator(updated_geometry2, filename2)
                        
                        gsub.write(f'sgx16 '+filename1+' 119 \n')
                        gsub.write(f'sgx16 '+filename2+' 119 \n')
                        gsub.write('\n')

                os.system(f"sbatch "+reduced_filename+".sub")
    print("Launch complete...") 
                    
    
logfilelist=[]
errorfiles=[]
for file in os.listdir():
    if "IRC" in file and "log" in file and not "logfile" in file:
        with open(file, 'r') as f:
            lines = f.readlines()
            if "Normal termination" in lines[-1]:
                logfilelist.append(file)
            else:
                errorfiles.append(file)

if len(errorfiles)!=0:
    answer=input("The opt+freq of the CREST conformers have led to some errors, do you wish to continue (yes or no)? ")
    if answer.lower()=='no':
        print("Exiting the program...")
        sys.exit()

launcher(logfilelist)
