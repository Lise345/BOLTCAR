import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


#with open('./parameters.txt', 'r') as parameters:
    #file_content = parameters.read()

    #energy_of_separate_reagents_match = re.search(r'Energy of separate reagents (.+)', file_content)
    #energy_of_separate_reagents = float(energy_of_separate_reagents_match.group(1))


def read_parameters(file_path):
    with open(file_path, 'r') as file:
        file_content = file.read()
    basis_1 = re.search(r'Basis (.+)', file_content)
    if basis_1:
        basis_1 = basis_1.group(1).strip()
    else:
        basis_1 = None
    return basis_1



def extract_values(file_path):
    pvdz_energy, pvtz_energy, pvqz_energy, gibbs_free_energy = None, None, None, None
    last_scf_done = None

    found_pvdz, found_pvtz, found_pvqz = False, False, False  # Track section existence

    print(f"Processing file: {file_path}")

    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if 'Error termination via Lnk1e' in line:
                return np.nan, np.nan, np.nan, np.nan  # Mark calculation as failed
            
            if 'SCF Done:  E(RM062X) =' in line:
                match = re.search(r'SCF Done:  E\(RM062X\) =\s+(-?\d+\.\d+)', line)
                if match:
                    last_scf_done = float(match.group(1))  # Always keep last SCF value

                    # Assign SCF Done value only if we are in the correct section
                    if found_pvdz and pvdz_energy is None:
                        pvdz_energy = last_scf_done
                    if found_pvtz and pvtz_energy is None:
                        pvtz_energy = last_scf_done
                    if found_pvqz and pvqz_energy is None:
                        pvqz_energy = last_scf_done

            # Detect basis set sections
            if 'm062x cc-pvdz' in line:
                found_pvdz = True
            if 'm062x cc-pvtz' in line:
                found_pvtz = True
            if 'm062x cc-pvqz' in line:
                found_pvqz = True

            if 'Thermal correction to Gibbs Free Energy=' in line:
                try:
                    gibbs_free_energy = float(line.split()[-1])
                except ValueError:
                    gibbs_free_energy = np.nan

    # Ensure missing values are explicitly marked as NaN if a section was never found
    if not found_pvdz:
        pvdz_energy = None
    if not found_pvtz:
        pvtz_energy = None
    if not found_pvqz:
        pvqz_energy = None
    if gibbs_free_energy is None:
        gibbs_free_energy = None

    print(f"Final extracted values for {file_path}:")
    print(f"PVDZ: {pvdz_energy}, PVTZ: {pvtz_energy}, PVQZ: {pvqz_energy}, Gibbs Free Energy: {gibbs_free_energy}")

    return pvdz_energy, pvtz_energy, pvqz_energy, gibbs_free_energy


def jobid(filename):
    match = re.findall(r"(\d+)", filename) #Finds all sequences of numbers
    if match:
        for group in match:
            if len(group)==4: #if sequence of numbers is 4 units long then return that as the jobid; unlikely that in name a sequence of 4 numbers would be added
                numbers=group
    return numbers

# Read basis_1 from parameters.txt
parameters_file = 'parameters.txt'
if not os.path.exists(parameters_file):
    raise FileNotFoundError(f"The file '{parameters_file}' does not exist.")
basis_1 = read_parameters(parameters_file)
use_cbs_logic = basis_1.lower() == "cbs"

directory = './'
data_dict = {}

for filename in os.listdir(directory):
    if filename.endswith('Complex.log') or filename.endswith('Complex_R1.log') or filename.endswith('Complex_R2.log') or filename.endswith('Product.log') or filename.endswith('SP.log'):
        identification_number = filename.split('-')[1].split('_')[0]
        file_path = os.path.join(directory, filename)
        
        pvdz_energy, pvtz_energy, pvqz_energy, gibbs_free_energy = extract_values(file_path)
        
        print(f"{identification_number} gives")
        print(f"{pvdz_energy} {pvtz_energy} {pvqz_energy}")
        
        if pvdz_energy is None and pvtz_energy is None and pvqz_energy is None and gibbs_free_energy is None:
            continue  # Skip files that contain the error message
        
        if identification_number not in data_dict:
            data_dict[identification_number] = {'ID Number': identification_number}

        if use_cbs_logic:
            # Use original logic for CBS
            if filename.endswith('Complex.log'):
                data_dict[identification_number].update({
                    'Complex PVDZ Energy': pvdz_energy,
                    'Complex PVTZ Energy': pvtz_energy,
                    'Complex PVQZ Energy': pvqz_energy,
                    'Complex Gibbs Correction': gibbs_free_energy,
                })
            elif filename.endswith('Complex_R1.log'):
                data_dict[identification_number].update({
                    'ComplexR1 PVDZ Energy': pvdz_energy,
                    'ComplexR1 PVTZ Energy': pvtz_energy,
                    'ComplexR1 PVQZ Energy': pvqz_energy,
                    'ComplexR1 Gibbs Correction': gibbs_free_energy,
                })
            elif filename.endswith('Complex_R2.log'):
                data_dict[identification_number].update({
                    'ComplexR2 PVDZ Energy': pvdz_energy,
                    'ComplexR2 PVTZ Energy': pvtz_energy,
                    'ComplexR2 PVQZ Energy': pvqz_energy,
                    'ComplexR2 Gibbs Correction': gibbs_free_energy,
                })
            elif filename.endswith('SP.log'):
                data_dict[identification_number].update({
                    'TS PVDZ Energy': pvdz_energy,
                    'TS PVTZ Energy': pvtz_energy,
                    'TS PVQZ Energy': pvqz_energy,
                    'TS Gibbs Correction': gibbs_free_energy,
                })
            elif filename.endswith('Product.log'):
                data_dict[identification_number].update({
                    'Product PVDZ Energy': pvdz_energy,
                    'Product PVTZ Energy': pvtz_energy,
                    'Product PVQZ Energy': pvqz_energy,
                    'Product Gibbs Correction': gibbs_free_energy,
                })
        else:
            # New logic when basis_1 is not CBS
            complex_basis_energy_label = f'Complex {basis_1} Energy'
            ts_basis_energy_label = f'TS {basis_1} Energy'
            product_basis_energy_label = f'Product {basis_1} Energy'

            if filename.endswith('Complex.log'):
                data_dict[identification_number].update({
                    complex_basis_energy_label: pvqz_energy,  # Assign PVQZ energy to Complex basis_1 energy
                    'Complex Gibbs Correction': gibbs_free_energy,
                })
            elif filename.endswith('SP.log'):
                data_dict[identification_number].update({
                    ts_basis_energy_label: pvqz_energy,  # Assign PVQZ energy to TS basis_1 energy
                    'TS Gibbs Correction': gibbs_free_energy,
                })
            elif filename.endswith('Product.log'):
                data_dict[identification_number].update({
                    product_basis_energy_label: pvqz_energy,  # Assign PVQZ energy to Product basis_1 energy
                    'Product Gibbs Correction': gibbs_free_energy,
                })

# Create DataFrame dynamically
df = pd.DataFrame.from_dict(data_dict, orient='index')




if use_cbs_logic:
    df['Extrapolated Complex Energy'] = (
        (df['Complex PVDZ Energy'] * df['Complex PVQZ Energy'] - df['Complex PVTZ Energy']**2)
        / (df['Complex PVDZ Energy'] + df['Complex PVQZ Energy'] - 2 * df['Complex PVTZ Energy'])
    )
    df['Extrapolated Reagent 1 Energy'] = (
        (df['ComplexR1 PVDZ Energy'] * df['ComplexR1 PVQZ Energy'] - df['ComplexR1 PVTZ Energy']**2)
        / (df['ComplexR1 PVDZ Energy'] + df['ComplexR1 PVQZ Energy'] - 2 * df['ComplexR1 PVTZ Energy'])
    )
    df['Extrapolated Reagent 2 Energy'] = (
        (df['ComplexR2 PVDZ Energy'] * df['ComplexR2 PVQZ Energy'] - df['ComplexR2 PVTZ Energy']**2)
        / (df['ComplexR2 PVDZ Energy'] + df['ComplexR2 PVQZ Energy'] - 2 * df['ComplexR2 PVTZ Energy'])
    )
    df['Extrapolated TS Energy'] = (
        (df['TS PVDZ Energy'] * df['TS PVQZ Energy'] - df['TS PVTZ Energy']**2)
        / (df['TS PVDZ Energy'] + df['TS PVQZ Energy'] - 2 * df['TS PVTZ Energy'])
    )
    df['Extrapolated Product Energy'] = (
        (df['Product PVDZ Energy'] * df['Product PVQZ Energy'] - df['Product PVTZ Energy']**2)
        / (df['Product PVDZ Energy'] + df['Product PVQZ Energy'] - 2 * df['Product PVTZ Energy'])
    )
else:
    complex_basis_energy_label = f'Complex {basis_1} Energy'
    ts_basis_energy_label = f'TS {basis_1} Energy'
    product_basis_energy_label = f'Product {basis_1} Energy'

    df['Extrapolated Complex Energy'] = df[complex_basis_energy_label]
    df['Extrapolated TS Energy'] = df[ts_basis_energy_label]
    df['Extrapolated Product Energy'] = df[product_basis_energy_label]

df['energy_of_separate_reagents'] = df['Extrapolated Reagent 1 Energy']+df['Extrapolated Reagent 2 Energy']+df['ComplexR1 Gibbs Correction']+df['ComplexR2 Gibbs Correction']
df['Minimal energy_of_separate_reagents'] = df['energy_of_separate_reagents'].min()

df['Complex Energy'] = 627.5 * (df['Extrapolated Complex Energy'] + df['Complex Gibbs Correction'] - df['Minimal energy_of_separate_reagents'])
df['TS Energy'] = 627.5 * (df['Extrapolated TS Energy'] + df['TS Gibbs Correction'] - df['Minimal energy_of_separate_reagents'])
df['Product Energy'] = 627.5 * (df['Extrapolated Product Energy'] + df['Product Gibbs Correction'] - df['Minimal energy_of_separate_reagents'])

# Ensure 'Complex Energy' is numeric and handle NaN or non-numeric values
df['Complex Energy'] = pd.to_numeric(df['Complex Energy'], errors='coerce')
df['TS Energy'] = pd.to_numeric(df['TS Energy'], errors='coerce').fillna(0)
df['Product Energy'] = pd.to_numeric(df['Product Energy'], errors='coerce').fillna(0)


# List of required energy columns for Pi calculation
required_columns = ['Complex PVDZ Energy', 'Complex PVTZ Energy', 'Complex PVQZ Energy', 'TS PVDZ Energy', 'TS PVTZ Energy', 'TS PVQZ Energy', 'Product PVDZ Energy', 'Product PVTZ Energy', 'Product PVQZ Energy']

# Set Complex Energy to NaN if any required energy value is missing
df.loc[df[required_columns].isnull().any(axis=1), ['Complex Energy', 'TS Energy', 'Product Energy']] = np.nan

# Compute Pi Value numerically where all values exist
df['Pi Value'] = np.exp(-df['Complex Energy'] / (0.001987204259 * 298.15))

# Ensure Pi Value is zero if Complex Energy > 1
df.loc[df['Complex Energy'] > 1, 'Pi Value'] = 0

# Create a separate column for display in Excel
df['Pi Value Display'] = df['Pi Value']

# Replace NaN values with "Calculation failed" ONLY in the Excel output column
df['Pi Value Display'] = df['Pi Value Display'].apply(lambda x: 'Calculation failed' if pd.isna(x) else x)





# Handle Pi Value sums safely
pi_sum = df['Pi Value'].sum(skipna=True) 

if pi_sum == 0:  # Avoid division by zero
    df['Percentage'] = 0
else:
    df['Percentage'] = df['Pi Value'].apply(lambda x: x / pi_sum if pd.notna(x) else 0)


# Calculate Rate Constant
df['Rate Constant'] = ((298.15 * 1.380649E-23) / 6.62607015E-34) * np.exp(-(df['TS Energy']) * 1000 * 4.184 / (8.314 * 298.15))

# Sort by identification number
df = df.sort_values(by='ID Number')

# Check for stable complexes
if (df['Complex Energy'] > 1).all():
    min_ts_energy_row = df.loc[df['TS Energy'].idxmin()]
    print(f"No stable complexes found, we assume the reaction will be the one with the lowest barrier, being {min_ts_energy_row['ID Number']}.")

# Save results
df.to_excel('FASTCAR_results.xlsx', index=False, columns=[col if col != 'Pi Value' else 'Pi Value Display' for col in df.columns])


print("Data extraction complete. The results are saved in 'FASTCAR_results.xlsx'.")

# Check for stable complexes
if (df['Complex Energy'] > 1).all():
    min_ts_energy_row = df.loc[df['TS Energy'].idxmin()]
    print(f"No stable Complexes found, the path of the lowest transition state energy will be followed, being {min_ts_energy_row['ID Number']}.")




# Convert Pi Value to numeric where possible (ignoring "Calculation failed")
df['Pi Value Numeric'] = pd.to_numeric(df['Pi Value'], errors='coerce')

# Check if there are any positive Pi Values
has_positive_pi = (df['Pi Value Numeric'] > 0).any()

# If there are positive Pi Values, filter for them; otherwise, use all IDs
if has_positive_pi:
    df_filtered = df[df['Pi Value Numeric'] > 0]
else:
    df_filtered = df  # Use all IDs if no positive Pi Values exist

# Ensure valid numeric values before plotting
df_filtered = df_filtered.dropna(subset=['Complex Energy', 'TS Energy', 'Product Energy'])




# Ensure the directory for saving plots exists
output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)

# Plot and save Complex Energy
plt.figure(figsize=(10, 5))
plt.scatter(df_filtered['ID Number'], df_filtered['Complex Energy'], color='b', label='Complex Energy')
plt.xlabel('ID Number')
plt.ylabel('Energy (kcal/mol)')
plt.title('Complex Energies for ID Numbers with Pi > 0')
plt.xticks(rotation=90)
plt.legend()
plt.grid()
plt.savefig(os.path.join(output_dir, "complex_energies.png"), dpi=300, bbox_inches='tight')
plt.close()

# Plot and save TS Energy
plt.figure(figsize=(10, 5))
plt.scatter(df_filtered['ID Number'], df_filtered['TS Energy'], color='r', label='TS Energy')
plt.xlabel('ID Number')
plt.ylabel('Energy (kcal/mol)')
plt.title('Transition State Energies for ID Numbers with Pi > 0')
plt.xticks(rotation=90)
plt.legend()
plt.grid()
plt.savefig(os.path.join(output_dir, "ts_energies.png"), dpi=300, bbox_inches='tight')
plt.close()

# Plot and save Product Energy
plt.figure(figsize=(10, 5))
plt.scatter(df_filtered['ID Number'], df_filtered['Product Energy'], color='g', label='Product Energy')
plt.xlabel('ID Number')
plt.ylabel('Energy (kcal/mol)')
plt.title('Product Energies for ID Numbers with Pi > 0')
plt.xticks(rotation=90)
plt.legend()
plt.grid()
plt.savefig(os.path.join(output_dir, "product_energies.png"), dpi=300, bbox_inches='tight')
plt.close()

# Plot and save Percentages
plt.figure(figsize=(10, 5))
plt.scatter(df_filtered['ID Number'], df_filtered['Percentage']*100, color='y', label='Percentages')
plt.xlabel('ID Number')
plt.ylabel('Percentage (%)')
plt.title('Percentages for ID Numbers with Pi > 0')
plt.xticks(rotation=90)
plt.legend()
plt.grid()
plt.savefig(os.path.join(output_dir, "percentages.png"), dpi=300, bbox_inches='tight')
plt.close()

print(f"Plots saved in '{output_dir}' directory.")

