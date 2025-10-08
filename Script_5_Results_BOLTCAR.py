import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


#with open('./parameters.txt', 'r') as parameters:
    #file_content = parameters.read()

    #structurename = re.search(r'Energy of separate reagents (.+)', file_content)
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
            
            if 'SCF Done:' in line:
                match = re.search(r'SCF Done:\s+E\(\S+\)\s+=\s+(-?\d+\.\d+)', line)
                if match:
                    last_scf_done = float(match.group(1))  # Always keep last SCF value

                    # Assign SCF Done value only if we are in the correct section
                    if found_pvdz and pvdz_energy is None:
                        pvdz_energy = last_scf_done
                    if found_pvtz and pvtz_energy is None:
                        pvtz_energy = last_scf_done
                    if found_pvqz:
                        pvqz_energy = last_scf_done

            # Detect basis set sections
            if 'opt=calcfc' and 'pvdz' in line:
                found_pvdz = True
            if 'opt=calcfc' and 'pvtz' in line:
                found_pvtz = True
            if 'opt=calcfc' and 'pvqz' in line:
                found_pvqz = True
            else:
                if 'opt=calcfc':
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
    if filename.endswith('Complex.log') or filename.endswith('Complex_R1.log') or filename.endswith('Complex_R2.log') or filename.endswith('Product.log') or filename.endswith('SP.log') or ('IRC' not in filename and filename.endswith('.log')):
        base_name = os.path.splitext(filename)[0]  # Removes .log
        parts = base_name.split('-')

        if len(parts) > 1:
            identification_number = parts[1].split('_')[0]
        else:
            print(f"⚠️ Could not extract ID from filename: {filename}")
            continue
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
            complex_basis_energy_label = f'Complex PVQZ Energy'
            ts_basis_energy_label = f'TS PVQZ Energy'
            product_basis_energy_label = f'Product PVQZ Energy'

            if filename.endswith('Complex.log'):
                data_dict[identification_number].update({
                    complex_basis_energy_label: pvqz_energy,  # Assign PVQZ energy to Complex basis_1 energy
                    'Complex Gibbs Correction': gibbs_free_energy,
                })
            elif 'IRC' not in filename and filename.endswith('.log'):
                data_dict[identification_number].update({
                    ts_basis_energy_label: pvqz_energy,  # Assign PVQZ energy to TS basis_1 energy
                    'TS Gibbs Correction': gibbs_free_energy,
                })
            elif filename.endswith('Product.log'):
                data_dict[identification_number].update({
                    product_basis_energy_label: pvqz_energy,  # Assign PVQZ energy to Product basis_1 energy
                    'Product Gibbs Correction': gibbs_free_energy,
                })
            elif filename.endswith('Complex_R1.log'):
                data_dict[identification_number].update({
                    'ComplexR1 PVQZ Energy': pvqz_energy,  # Assign PVQZ energy to ComplexR1 basis_1 energy
                    'ComplexR1 Gibbs Correction': gibbs_free_energy,
                })
            elif filename.endswith('Complex_R2.log'):
                data_dict[identification_number].update({
                    'ComplexR2 PVQZ Energy': pvqz_energy,  # Assign PVQZ energy to ComplexR2 basis_1 energy
                    'ComplexR2 Gibbs Correction': gibbs_free_energy,
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
    complex_basis_energy_label = f'Complex PVQZ Energy'
    ts_basis_energy_label = f'TS PVQZ Energy'
    product_basis_energy_label = f'Product PVQZ Energy'

    df['Extrapolated Complex Energy'] = df[complex_basis_energy_label]
    df['Extrapolated TS Energy'] = df[ts_basis_energy_label]
    df['Extrapolated Product Energy'] = df[product_basis_energy_label]
    df['Extrapolated Reagent 1 Energy'] = df.get('ComplexR1 PVQZ Energy', np.nan)
    df['Extrapolated Reagent 2 Energy'] = df.get('ComplexR2 PVQZ Energy', np.nan)

df['energy_of_separate_reagents'] = df['Extrapolated Reagent 1 Energy']+df['Extrapolated Reagent 2 Energy']+df['ComplexR1 Gibbs Correction']+df['ComplexR2 Gibbs Correction']
df['Minimal energy_of_separate_reagents'] = df['energy_of_separate_reagents'].min()
df['delta energy of separate reagents'] = df['energy_of_separate_reagents'] - df['Minimal energy_of_separate_reagents']

df['Complex Energy'] = 627.5 * (df['Extrapolated Complex Energy'] + df['Complex Gibbs Correction'] - df['Minimal energy_of_separate_reagents'])
df['TS Energy'] = 627.5 * (df['Extrapolated TS Energy'] + df['TS Gibbs Correction'] - df['Minimal energy_of_separate_reagents'])
df['Product Energy'] = 627.5 * (df['Extrapolated Product Energy'] + df['Product Gibbs Correction'] - df['Minimal energy_of_separate_reagents'])

# Ensure 'Complex Energy' is numeric and handle NaN or non-numeric values
df['Complex Energy'] = pd.to_numeric(df['Complex Energy'], errors='coerce')
df['TS Energy'] = pd.to_numeric(df['TS Energy'], errors='coerce').fillna(0)
df['Product Energy'] = pd.to_numeric(df['Product Energy'], errors='coerce').fillna(0)


if use_cbs_logic:
    # List of required energy columns for Pi calculation
    required_columns = ['Complex PVDZ Energy', 'Complex PVTZ Energy', 'Complex PVQZ Energy', 'TS PVDZ Energy', 'TS PVTZ Energy', 'TS PVQZ Energy', 'Product PVDZ Energy', 'Product PVTZ Energy', 'Product PVQZ Energy']
else:
    # List of required energy columns for Pi calculation when not using CBS logic
    required_columns = ['Extrapolated Complex Energy', 'Extrapolated TS Energy', 'Extrapolated Product Energy']

# Set Complex Energy to NaN if any required energy value is missing
df.loc[df[required_columns].isnull().any(axis=1), ['Complex Energy', 'TS Energy', 'Product Energy']] = np.nan


### Calculating percentages for Boltzmann ###

# Constants
R = 0.001987204259   # kcal·mol^-1·K^-1
T = 298.15
RT = R * T


# Ensure numeric
df['Complex Energy'] = pd.to_numeric(df['Complex Energy'], errors='coerce')
df['delta energy of separate reagents'] = pd.to_numeric(df['delta energy of separate reagents'], errors='coerce')

# Compute Pi Value forward per requirement:
neg_mask = df['Complex Energy'] < 0

if neg_mask.any():
    # If at least one Complex Energy < 0:
    #   - rows with Complex Energy < 0:  exp(-ComplexEnergy/RT)
    #   - rows with Complex Energy >= 0 or NaN: 0
    df['Pi Value forward'] = 0.0
    df.loc[neg_mask, 'Pi Value forward'] = np.exp(-df.loc[neg_mask, 'Complex Energy'] / RT)
else:
    # Otherwise fall back to the original definition
    df['Pi Value forward'] = np.exp(-df['delta energy of separate reagents'] / RT)
    
df['Pi Value reverse'] = np.exp(-df['Product Energy'] / (0.001987204259 * 298.15))

# Create a separate column for display in Excel
df['Pi Value Display F'] = df['Pi Value forward']
df['Pi Value Display R'] = df['Pi Value reverse']

# Replace NaN values with "Calculation failed" ONLY in the Excel output column
df['Pi Forward'] = df['Pi Value Display F'].apply(lambda x: 'Calculation failed' if pd.isna(x) else x)
df['Pi Reverse'] = df['Pi Value Display R'].apply(lambda x: 'Calculation failed' if pd.isna(x) else x)

# Convert Pi columns to numeric safely (invalid entries become NaN)
df['Pi Forward Numeric'] = pd.to_numeric(df['Pi Forward'], errors='coerce')
df['Pi Reverse Numeric'] = pd.to_numeric(df['Pi Reverse'], errors='coerce')

# Recalculate sums with only valid numeric entries
pi_sum_forward = df['Pi Forward Numeric'].sum(skipna=True)
pi_sum_reverse = df['Pi Reverse Numeric'].sum(skipna=True)

if pi_sum_forward == 0:  # Avoid division by zero
    df['Percentage Forward'] = 0
else:
    df['Percentage Forward'] = df['Pi Forward Numeric'].apply(lambda x: x / pi_sum_forward if pd.notna(x) else 0).round(3)*100

if pi_sum_reverse == 0:
    df['Percentage Reverse'] = 0
else:
    df['Percentage Reverse'] = df['Pi Reverse Numeric'].apply(lambda x: x / pi_sum_reverse if pd.notna(x) else 0).round(3)*100
    

df['Forward Barrier'] = df['TS Energy'] - df['delta energy of separate reagents']
df['Reverse Barrier'] = df['TS Energy'] - df['Product Energy']


# Calculate Rate Constant
df['Forward Rate Constant'] = ((298.15 * 1.380649E-23) / 6.62607015E-34) * np.exp(-(df['Forward Barrier']) * 1000 * 4.184 / (8.314 * 298.15))
df['Reverse Rate Constant'] = ((298.15 * 1.380649E-23) / 6.62607015E-34) * np.exp(-(df['Reverse Barrier']) * 1000 * 4.184 / (8.314 * 298.15))

avg_forward = "Average not relevant, fastest reaction shown below"
avg_reverse = "Error: no stable adducts found"

# Ensure valid weights and rates for forward
valid_forward_mask = pd.notna(df['Forward Rate Constant']) & pd.notna(df['Percentage Forward']) & (df['Percentage Forward'] > 0)
if valid_forward_mask.any():
    avg_forward_value = np.average(df.loc[valid_forward_mask, 'Forward Rate Constant'],
                                   weights=df.loc[valid_forward_mask, 'Percentage Forward'])
    avg_forward = "{:.1E}".format(avg_forward_value)

# Ensure valid weights and rates for reverse
valid_reverse_mask = pd.notna(df['Reverse Rate Constant']) & pd.notna(df['Percentage Reverse']) & (df['Percentage Reverse'] > 0)
if valid_reverse_mask.any():
    avg_reverse_value = np.average(df.loc[valid_reverse_mask, 'Reverse Rate Constant'],
                                   weights=df.loc[valid_reverse_mask, 'Percentage Reverse'])
    avg_reverse = "{:.1E}".format(avg_reverse_value)



lowest_forward = ""

# Check if the avg_forward message was used
if avg_forward == "Average not relevant, fastest reaction shown in Excel":
    # Instead of relying on valid_forward_mask, just use all valid Forward Rate Constants
    valid_forward_rate_mask = pd.notna(df['Forward Rate Constant']) & (df['Forward Rate Constant'] > 0)
    if valid_forward_rate_mask.any():
        lowest_forward_value = df.loc[valid_forward_rate_mask, 'Forward Rate Constant'].max()
        lowest_forward = "{:.1E}".format(lowest_forward_value)


avg_data = pd.DataFrame({
    'Weighted Average Forward rate Constant': [avg_forward, lowest_forward if lowest_forward else None],
    'Weighted Average Reverse rate Constant': [avg_reverse, None]
})

# Sort by identification number
df = df.sort_values(by='ID Number')

# Check for stable complexes
if (df['Complex Energy'] > 1).all():
    min_ts_energy_row = df.loc[df['TS Energy'].idxmin()]
    print(f"No stable complexes found, we assume the reaction will be the one with the lowest barrier, being {min_ts_energy_row['ID Number']}.")

round_cols = ['Complex Energy', 'TS Energy', 'Product Energy', 'Pi Value Display', 'Percentage', 'Forward Barrier', 'Reverse Barrier']
for col in round_cols:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce').round(1)


first_cols = ['ID Number', 'Complex Energy', 'TS Energy', 'Product Energy', 'Forward Barrier', 'Percentage Forward', 'Reverse Barrier', 'Percentage Reverse', 'Forward Rate Constant', 'Reverse Rate Constant']

# Then add the remaining ones that are not already in first_cols
remaining_cols = [col for col in df.columns if col not in first_cols]

# Full column order
new_columns = first_cols + remaining_cols

# Save results
with pd.ExcelWriter('BOLTCAR_results.xlsx') as writer:
    # Save the full dataset
    df.to_excel(writer, index=False, sheet_name='Full Results', columns=[col if col != 'Pi Value' else 'Pi Value Display' for col in new_columns])

    # Save the weighted averages
    avg_data.to_excel(writer, index=False, sheet_name='Weighted Averages')

print("Data extraction complete. The results are saved in 'BOLTCAR_results.xlsx'.")


# Convert Pi Value to numeric where possible (ignoring "Calculation failed")
df['Pi Value Numeric'] = pd.to_numeric(df['Pi Value forward'], errors='coerce')

# Check if there are any positive Pi Values
has_positive_pi = (df['Pi Value Numeric'] > 0).any()

# If there are positive Pi Values, filter for them; otherwise, use all IDs
if has_positive_pi:
    df_filtered = df[df['Pi Value Numeric'] > 0]
else:
    df_filtered = df  # Use all IDs if no positive Pi Values exist

# Ensure valid numeric values before plotting
df_filtered = df_filtered.dropna(subset=['Complex Energy', 'TS Energy', 'Product Energy'])

# Define y-limits with margin
min_energy = min(df_filtered['Complex Energy'].min(), df_filtered['TS Energy'].min(), df_filtered['Product Energy'].min()) - 1
max_energy = max(df_filtered['Complex Energy'].max(), df_filtered['TS Energy'].max(), df_filtered['Product Energy'].max()) + 1

# Ensure the directory for saving plots exists
output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)

# Plot and save Complex Energy
plt.figure(figsize=(10, 5))
plt.scatter(df_filtered['ID Number'], df_filtered['Complex Energy'], color='b', label=f'Complex Energy\nAvg kf = {avg_forward}')
plt.xlabel('ID Number')
plt.ylabel('Energy (kcal/mol)')
plt.title('Complex Energies for ID Numbers with Pi > 0')
plt.xticks(rotation=90)
plt.legend()
plt.ylim(min_energy, max_energy)
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
plt.ylim(min_energy, max_energy)
plt.grid()
plt.savefig(os.path.join(output_dir, "ts_energies.png"), dpi=300, bbox_inches='tight')
plt.close()

# Plot and save Product Energy
plt.figure(figsize=(10, 5))
plt.scatter(df_filtered['ID Number'], df_filtered['Product Energy'], color='g', label=f'Product Energy\nAvg kr = {avg_reverse}')
plt.xlabel('ID Number')
plt.ylabel('Energy (kcal/mol)')
plt.title('Product Energies for ID Numbers with Pi > 0')
plt.xticks(rotation=90)
plt.legend()
plt.ylim(min_energy, max_energy)
plt.grid()
plt.savefig(os.path.join(output_dir, "product_energies.png"), dpi=300, bbox_inches='tight')
plt.close()

# Plot and save Percentages
plt.figure(figsize=(10, 5))
plt.scatter(df_filtered['ID Number'], df_filtered['Percentage Forward'], color='orange', label='Forward %')
plt.scatter(df_filtered['ID Number'], df_filtered['Percentage Reverse'], color='purple', label='Reverse %')
plt.xlabel('ID Number')
plt.ylabel('Percentage (%)')
plt.title('Percentages for ID Numbers with Pi > 0')
plt.xticks(rotation=90)
plt.legend()
plt.grid()
plt.savefig(os.path.join(output_dir, "percentages.png"), dpi=300, bbox_inches='tight')
plt.close()




# Create 2x2 subplot layout
fig, axs = plt.subplots(2, 2, figsize=(16, 12))

# Plot 1: Complex Energy
axs[0, 0].scatter(df_filtered['ID Number'], df_filtered['Complex Energy'], color='b', label=f'Complex Energy\nAvg kf = {avg_forward}')
axs[0, 0].set_title('Complex Energies')
axs[0, 0].set_xlabel('ID Number')
axs[0, 0].set_ylabel('Energy (kcal/mol)')
axs[0, 0].set_ylim(min_energy, max_energy)
axs[0, 0].legend()
axs[0, 0].grid(True, which='both', linestyle='--', color='lightgrey')
axs[0, 0].tick_params(axis='x', rotation=90)

# Plot 2: TS Energy
axs[0, 1].scatter(df_filtered['ID Number'], df_filtered['TS Energy'], color='r', label='TS Energy')
axs[0, 1].set_title('Transition State Energies')
axs[0, 1].set_xlabel('ID Number')
axs[0, 1].set_ylabel('Energy (kcal/mol)')
axs[0, 1].set_ylim(min_energy, max_energy)
axs[0, 1].legend()
axs[0, 1].grid(True, which='both', linestyle='--', color='lightgrey')
axs[0, 1].tick_params(axis='x', rotation=90)

# Plot 3: Product Energy
axs[1, 0].scatter(df_filtered['ID Number'], df_filtered['Product Energy'], color='g', label=f'Product Energy\nAvg kr = {avg_reverse}')
axs[1, 0].set_title('Product Energies')
axs[1, 0].set_xlabel('ID Number')
axs[1, 0].set_ylabel('Energy (kcal/mol)')
axs[1, 0].set_ylim(min_energy, max_energy)
axs[1, 0].legend()
axs[1, 0].grid(True, which='both', linestyle='--', color='lightgrey')
axs[1, 0].tick_params(axis='x', rotation=90)

# Plot 4: Percentages
axs[1, 1].scatter(df_filtered['ID Number'], df_filtered['Percentage Forward'], color='orange', label='Forward %')
axs[1, 1].scatter(df_filtered['ID Number'], df_filtered['Percentage Reverse'], color='purple', label='Reverse %')
axs[1, 1].set_title('Reaction Percentages')
axs[1, 1].set_xlabel('ID Number')
axs[1, 1].set_ylabel('Percentage (%)')
axs[1, 1].legend()
axs[1, 1].grid(True, which='both', linestyle='--', color='lightgrey')
axs[1, 1].tick_params(axis='x', rotation=90)

plt.tight_layout()

# Save as PDF
pdf_path = os.path.join(output_dir,"BOLTCAR_combined_plots.pdf")
plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"Combined 2x2 plots saved as '{pdf_path}'.")

print(f"Plots saved in '{output_dir}' directory.")
