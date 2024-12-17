import os
import re

# Atomic symbols mapping
atomic_symbols = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca",
    21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn",
    31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr"
}

def extract_coordinates_from_log(file_path):
    """Extracts atomic symbols and coordinates from the LAST Standard orientation block."""
    with open(file_path, 'r') as file:
        lines = file.readlines()

    block_starts = [i + 5 for i, line in enumerate(lines) if "Standard orientation" in line]

    if not block_starts:
        print(f"Warning: 'Standard orientation' block not found in {file_path}")
        return []

    start = block_starts[-1]  # Take the last occurrence
    end = start

    for i in range(start, len(lines)):
        if "-----" in lines[i]:
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

def write_gaussian_input(file_name, molecule, suffix):
    """Writes a Gaussian input file with the extracted coordinates."""
    header = f"""%nprocshared=8
%mem=16GB
%chk={file_name}.chk
# opt=calcfc freq m062x cc-pvdz empiricaldispersion=gd3 

{file_name}_{suffix} optfreq

0 1\n"""

    output_path = f"{file_name}_{suffix}.gjf"
    with open(output_path, 'w') as output_file:
        output_file.write(header)
        for atom in molecule:
            output_file.write(f" {atom[0]:<2} {atom[1]:>15.8f} {atom[2]:>15.8f} {atom[3]:>15.8f}\n")
        output_file.write("\n")
    print(f"Gaussian input file written: {output_path}")

def process_all_logs(directory, parameters_file):
    """Processes all *Complex.log files in the directory."""
    molecule1_indices, molecule2_indices = read_parameters(parameters_file)

    for filename in os.listdir(directory):
        if filename.endswith("Complex.log"):
            file_path = os.path.join(directory, filename)
            print(f"Processing file: {filename}")

            extracted_atoms = extract_coordinates_from_log(file_path)
            if not extracted_atoms:
                print(f"Skipping {filename}: No coordinates found.")
                continue

            # Ensure indices are valid
            max_index = len(extracted_atoms)
            if any(i > max_index for i in molecule1_indices + molecule2_indices):
                print(f"Error: Indices exceed available atoms in {filename}.")
                continue

            # Extract molecule coordinates
            molecule1 = [extracted_atoms[i - 1] for i in molecule1_indices]
            molecule2 = [extracted_atoms[i - 1] for i in molecule2_indices]

            # Generate Gaussian input files for R1 and R2
            base_name = os.path.splitext(filename)[0]
            write_gaussian_input(base_name, molecule1, "R1")
            write_gaussian_input(base_name, molecule2, "R2")

def main():
    directory = "./"  # Directory with .log files
    parameters_file = "parameters.txt"  # Path to parameters.txt

    process_all_logs(directory, parameters_file)

if __name__ == "__main__":
    main()
