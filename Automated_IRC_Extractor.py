import os
import math
import sys
import re
import openpyxl

with open('./parameters.txt', 'r') as parameters:
    file_content = parameters.read()
    print("File content read by script:")
    print(repr(file_content))  # Use repr() to reveal hidden characters like \n or spaces

    time_statp = re.search(r'Time for stationary calcs\s+(\d+)', file_content)
    if time_statp:
        time_for_stationary_calcs = int(time_statp.group(1))
        statp_time = f'{time_for_stationary_calcs}:00:00'  # Format to HH:MM:SS
    else:
        statp_time = '25:00:00'  # Default value if not found
        print("Time for stationary calculations not found, defaulting to 25:00:00")

    time_TS = re.search(r'Time for TS calcs\s+(\d+)', file_content)
    if time_TS:
        time_for_TS = int(time_TS.group(1))
        TS_time = f'{time_for_TS}:00:00'  # Format to HH:MM:SS
    else:
        TS_time = '25:00:00'  # Default value if not found
        print("Time for stationary calculations not found, defaulting to 25:00:00")

    # Extract the size_molecule value
    size_molecule_match = re.search(r'size_molecule\s*=\s*(\d+)', file_content)
    if size_molecule_match:
        size_molecule = int(size_molecule_match.group(1))

    # Extract the CC1_in value
    CC1_in_match = re.search(r'CC1_in\s*=\s*"(.*?)"', file_content)
    if CC1_in_match:
        CC1_in = list(map(int, CC1_in_match.group(1).split()))

    rootdir_match = re.search(r'rootdir\s+(.+)', file_content)
    rootdir = rootdir_match.group(1).strip().strip("'\"")

    binfolder = re.search(r'bin (.+)', file_content)
    binfolder = binfolder.group(1)

    #DFT parameters
    basis_in= re.search(r'Basis (.+)', file_content)
    basis_in= basis_in.group(1).strip()
    if basis_in.lower()=='cbs':
        basis_1='cc-pvdz'
        basis_2='cc-pvtz'
        basis_3='cc-pvqz'
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
    if charge:
        charge = charge.group(1)
    else:
        raise ValueError("Charge value not found in parameters file.")

    functional = re.search(r'Functional (.+)', file_content)
    functional = functional.group(1).strip()

    multiplicity = re.search(r'Multiplicity (-?\d+)', file_content)
    multiplicity = multiplicity.group(1)

CC1_out=CC1_in
CC1=CC1_in

#-------------Prepare TS geometries--------------

workbook = openpyxl.load_workbook("TS_analysis.xlsx")
sheet = workbook.active


data = []
for column in sheet.iter_cols(values_only=True):
    data.append(list(column))

print(data)

listofTS=[]
for caseofTS in data:
    if data.index(caseofTS)==3:
        for i in caseofTS:
            if caseofTS.index(i)>0 and i!=None:
                listofTS.append(i)

def SP_inputgenerator(xyzfile,filename):
    with open(xyzfile, 'r') as file:
        lines = file.readlines()
    with open(filename, 'x') as ip:
        ip.writelines("%nprocshared=8\n")
        ip.writelines("%mem=16GB\n")
        ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
        ip.writelines(f"# opt=(calcfc,ts,noeigentest) freq {functional} {basis_1} {dispersion} {solvent}\n")
        ip.writelines("\n")
        Title=filename[:-4]+" "+"cc-pVTZ_SP"+"\n"
        ip.writelines(Title)
        ip.writelines("\n")
        ip.writelines(f"{charge} {multiplicity}\n")
        
        for atom in lines[2:]:
            ip.writelines(atom)
        ip.writelines("\n")
        
        #Writing Link1 part for cc-pVTZ
        ip.writelines("--Link1--\n")
        ip.writelines("%nprocshared=8\n")
        ip.writelines("%mem=16GB\n")
        ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
        ip.writelines(f"# {functional} {basis_2} {dispersion} {solvent} Geom=Checkpoint \n")
        ip.writelines("\n")
        Title=filename[:-4]+" "+"cc-pVQT_SP"+"\n"
        ip.writelines(Title)
        ip.writelines("\n")
        ip.writelines(f"{charge} {multiplicity}\n")
        ip.writelines("\n")
        

	    #Writing Link1 part for cc-pVQZ
        ip.writelines("--Link1--\n")
        ip.writelines("%nprocshared=8\n")
        ip.writelines("%mem=16GB\n")
        ip.writelines("%chk="+filename[:-4]+".chk"+"\n")
        ip.writelines(f"# {functional} {basis_3} {dispersion} {solvent} Geom=Checkpoint \n")
        ip.writelines("\n")
        Title=filename[:-4]+" "+"cc-pVQZ_SP"+"\n"
        ip.writelines(Title)
        ip.writelines("\n")
        ip.writelines(f"{charge} {multiplicity}\n")
        ip.writelines("\n")
        ip.close()
    return print("Input generated for " + filename[:-4])

#-------------Prepare Complex and Product geometries---------------

def geometryextractor(logfile):
    with open(logfile, 'r') as file:
        lines = file.readlines()
    #convergence_indices = []
    #for i, line in enumerate(lines):
    #    if "Delta-x Convergence Met" in line:
    #        convergence_indices.append(i)
    convergence_indices = []
    for i, line in enumerate(lines):
        if "Input orientation" in line:
            convergence_indices.append(i)

    conv_geom=convergence_indices[-1]+5
    
    
    #print("convergence found")
    
    #index=conv_geom
    #found=False
    #while found==False:
    #    if "Number     Number" in lines[index]:
    #        found=True
    #        index=index+2
    #    else:
    #        index=index-1
    #        continue
    #print("geometry found ...")
    
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
    #print(updated_geometry)
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
        
        #Writing Link1 part for cc-pVTZ
        if basis_in.lower()=='cbs':
	        ip.writelines("--Link1--\n")
	        ip.writelines("%nprocshared=8\n")
	        ip.writelines("%mem=16GB\n")
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
	        ip.writelines("%mem=16GB\n")
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


#-------------Launch calculations-----------------


import os
import re
import subprocess

inp_file_job_ids = []

def launcherstatp(logfilelist):
    for i, logfile1 in enumerate(logfilelist):
        number = logfile1[17:21]
        for j, logfile2 in enumerate(logfilelist):
            if number in logfile2 and logfile1 != logfile2 and j > i:
                geometry1 = geometryextractor(logfile1)
                updated_geometry1 = geometryconverter(geometry1)
                distance1_CC1 = distance(updated_geometry1, CC1)

                geometry2 = geometryextractor(logfile2)
                updated_geometry2 = geometryconverter(geometry2)
                distance2_CC1 = distance(updated_geometry2, CC1)

                reduced_filename = logfile1[:-4] + "_optE"
                with open(reduced_filename + ".sub", "w") as gsub:
                    gsub.write('#!/bin/sh\n')
                    gsub.write(f'#SBATCH --job-name={reduced_filename}\n')
                    gsub.write('#SBATCH --ntasks=12\n')
                    gsub.write(f'#SBATCH --output={reduced_filename}.logfile\n')
                    gsub.write(f'#SBATCH --time={statp_time}\n')
                    gsub.write('\n')
                    gsub.write('module load Gaussian/G16.A.03-intel-2022a\n')
                    gsub.write('export GAUSS_SCRDIR=$TMPDIR\n')
                    gsub.write('mkdir -p $GAUSS_SCRDIR\n')
                    
                    # Write Gaussian job commands
                    if distance1_CC1 > distance2_CC1:
                        filename1 = logfile1[:-4] + "_Complex.gjf"
                        filename2 = logfile2[:-4] + "_Product.gjf"
                        output_filename1 = filename1.replace(".gjf", ".log")
                        output_filename2 = filename2.replace(".gjf", ".log")
                    else:
                        filename2 = logfile2[:-4] + "_Complex.gjf"
                        filename1 = logfile1[:-4] + "_Product.gjf"
                        output_filename1 = filename1.replace(".gjf", ".log")
                        output_filename2 = filename2.replace(".gjf", ".log")                       

                    inputgenerator(updated_geometry1, filename1)
                    inputgenerator(updated_geometry2, filename2)
                    gsub.write(f'g16 < {filename1} > {output_filename1} &\n')
                    gsub.write(f'g16 < {filename2} > {output_filename2} &\n')
                    gsub.write(f'wait\n')
                    
                # Submit the job
                result = subprocess.run(
                    f"sbatch {reduced_filename}.sub",
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
                        print(f"Job submitted: {job_id}")
                else:
                    print(f"Failed to submit job: {result.stderr}")
          

def launcherTS(xyzlist):
    for xyzfile in xyzlist:
        filename=xyzfile[:-4]+"_SP.gjf"
        outputfilename=xyzfile[:-4]+"_SP.log"
        reduced_filename=filename[:-4]
        with open(reduced_filename+".sub","w") as gsub:
            gsub.write('#!/bin/sh\n')
            gsub.write(f'#SBATCH --job-name='+filename[:-4]+'\n')
            gsub.write('#SBATCH --ntasks=12\n')
            gsub.write(f'#SBATCH --output='+filename[:-4]+'.logfile\n')
            gsub.write(f'#SBATCH --time={TS_time}\n')
            gsub.write('\n')
            gsub.write('# Loading modules\n')
            gsub.write('module load Gaussian/G16.A.03-intel-2022a\n')  # Adjust based on the available Gaussian module
            gsub.write('\n')
            gsub.write('# Setting up Gaussian environment\n')
            gsub.write('export GAUSS_SCRDIR=$TMPDIR\n')  # Temporary directory for Gaussian scratch files
            gsub.write('mkdir -p $GAUSS_SCRDIR\n')
            gsub.write('#Launching calculation\n')
            gsub.write('export PATH={binfolder}:$PATH\n')
            
            SP_inputgenerator(xyzfile,filename)
            
            gsub.write(f'g16 < {filename} > {outputfilename}\n')
            gsub.write('\n')
            gsub.close()
        
# Submit the job
        result = subprocess.run(
            f"sbatch {reduced_filename}.sub",
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
                print(f"Job submitted: {job_id}")
        else:
            print(f"Failed to submit job: {result.stderr}") 



#-----dependent jobs -----
def launch_dependent_job():
    if inp_file_job_ids:
        dependency_str = ":".join(inp_file_job_ids)
        extractor_script = os.path.join(rootdir, '4_SeparateReagents.sub')
        
        if not os.path.exists(extractor_script):
            raise FileNotFoundError(f"Extractor script not found: {extractor_script}")

        dependency_command = [
            "sbatch",
            f"--dependency=afterany:{dependency_str}",
            extractor_script
        ]
        result = subprocess.run(
            dependency_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        if result.returncode == 0:
            print(f"Extractor job submitted successfully: {result.stdout}")
        else:
            print(f"Failed to submit extractor job: {result.stderr}")
    else:
        print("No jobs were submitted, skipping dependency job submission.")


#-----error logging-----
    
logfilelist=[]
errorfiles=[]
for file in os.listdir():
    if "IRC" in file and "log" in file and not "logfile" in file:
        with open(file, 'r') as f:
            lines = f.readlines()
            if " Pt  1 Step number   5 out of a maximum of" in lines[-1]:
                logfilelist.append(file)
            else:
                errorfiles.append(file)


launcherstatp(logfilelist)
if basis_in.lower()=='cbs':
	launcherTS(listofTS)
launch_dependent_job()

