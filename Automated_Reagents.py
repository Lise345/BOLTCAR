import os
import re
import numpy as np

# Atomic symbols mapping
atomic_symbols = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca",
    21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn",
    31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr"
}

def extract_coordinates_from_log(file_path):
    """Extracts atomic symbols and coordinates from a Complex.log file."""
    with open(file_path, 'r') as file:
        lines = file.readlines()

    start, end = None, None
    for i, line in enumerate(lines):
        if "Standard orientation" in line:
            start = i + 5  # Start after the dashed lines
        elif start and "-----" in line:
            end = i
            break

    extracted_data = []
    for line in lines[start:end]:
        parts = line.split()
        if len(parts) >= 6:
            atom_number = int(parts[1])
            x, y, z = map(float, parts[3:6])
            atom_symbol = atomic_symbols.get(atom_number, f"Unknown({atom_number})")
            extracted_data.append([atom_symbol, x, y, z])
    return extracted_data

def read_parameters(file_path):
    """Reads molecule atom indices from parameters.txt."""
    with open(file_path, 'r') as file:
        file_content = file.read()

    molecule1_atoms = list(map(int, re.search(r'molecule1_atoms\s*=\s*(.+)', file_content).group(1).split()))
    molecule2_atoms = list(map(int, re.search(r'molecule2_atoms\s*=\s*(.+)', file_content).group(1).split()))
    return molecule1_atoms, molecule2_atoms

def write_gaussian_input(file_name, molecule1, molecule2):
    """Writes a Gaussian input file with the extracted coordinates."""
    header = f"""%nprocshared=8
%mem=16GB
%chk={file_name}.chk
# opt=calcfc freq m062x cc-pvdz empiricaldispersion=gd3 

{file_name} optfreq

0 1\n"""

    output_path = f"{file_name}_R1.gjf"
    with open(output_path, 'w') as output_file:
        output_file.write(header)
        for atom in molecule1 + molecule2:
            output_file.write(f" {atom[0]:<2} {atom[1]:>15.8f} {atom[2]:>15.8f} {atom[3]:>15.8f}\n")

        output_file.write("\n--Link1--\n%nprocshared=8\n%mem=16GB\n%chk={file_name}.chk\n")
        output_file.write("# m062x cc-pvtz empiricaldispersion=gd3  Geom=Checkpoint\n\n")
        output_file.write(f"{file_name} optfreq E_ccpvtz\n\n0 1\n\n")
        output_file.write("--Link1--\n%nprocshared=8\n%mem=16GB\n%chk={file_name}.chk\n")
        output_file.write("# m062x cc-pvqz empiricaldispersion=gd3  Geom=Checkpoint\n\n")
        output_file.write(f"{file_name} optfreq E_ccpvqz\n\n0 1\n")

    print(f"Gaussian input file written: {output_path}")

def process_all_logs(directory, parameters_file):
    """Processes all *Complex.log files in the directory."""
    molecule1_indices, molecule2_indices = read_parameters(parameters_file)
    
    for filename in os.listdir(directory):
        if filename.endswith("Complex.log"):
            file_path = os.path.join(directory, filename)
            print(f"Processing file: {filename}")
            
            # Extract coordinates
            extracted_atoms = extract_coordinates_from_log(file_path)
            
            # Split into molecule 1 and molecule 2
            molecule1 = [extracted_atoms[i - 1] for i in molecule1_indices]
            molecule2 = [extracted_atoms[i - 1] for i in molecule2_indices]
            
            # Generate Gaussian input file
            base_name = os.path.splitext(filename)[0]
            write_gaussian_input(base_name, molecule1, molecule2)

def main():
    directory = "./"  # Set to the folder containing your .log files
    parameters_file = "parameters.txt"  # Path to parameters.txt
    
    process_all_logs(directory, parameters_file)

if __name__ == "__main__":
    main()
