import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# Constants
R = 8.314   # kcal·mol^-1·K^-1
T = 298.15
RT = R * T


#with open('./parameters.txt', 'r') as parameters:
    #file_content = parameters.read()

    #structurename = re.search(r'Energy of separate reagents (.+)', file_content)
    #Gibbs of separate reagents = float(Gibbs of separate reagents_match.group(1))


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
    import re
    import numpy as np
    from collections import deque

    # Regexes (allow D/E notation)
    scf_re = re.compile(r'SCF Done:\s+E\([^)]*\)\s*=\s*(-?\d+(?:\.\d+)?(?:[DEde][+\-]?\d+)?)')
    gibbs_re = re.compile(r'Thermal correction to Gibbs Free Energy\s*=\s*(-?\d+(?:\.\d+)?(?:[DEde][+\-]?\d+)?)')
    enth_re = re.compile(r' Thermal correction to Enthalpy\s*=\s*(-?\d+(?:\.\d+)?(?:[DEde][+\-]?\d+)?)')


    last_three_scf = deque(maxlen=3)  # will keep only the last 3 energies
    gibbs = None
    enthalpy = None

    with open(file_path, 'r', errors='ignore') as fh:
        for line in fh:
            low = line.lower()

            # Bail out on Gaussian fatal error
            if "error termination via lnk1e" in low:
                return np.nan, np.nan, np.nan, np.nan, np.nan

            m = scf_re.search(line)
            if m:
                try:
                    val = float(m.group(1).replace('D', 'E').replace('d', 'e'))
                    last_three_scf.append(val)
                except ValueError:
                    pass

            g = gibbs_re.search(line)
            if g:
                try:
                    gibbs = float(g.group(1).replace('D', 'E').replace('d', 'e'))
                except ValueError:
                    pass
            
            e = enth_re.search(line)
            if e:
                try:
                    enthalpy = float(e.group(1).replace('D', 'E').replace('d', 'e'))
                except ValueError:
                    pass

    # Map last three SCF energies -> (pVDZ, pVTZ, pVQZ)
    pvdz = pvtz = pvqz = np.nan
    n = len(last_three_scf)
    if n == 1:
        pvqz = last_three_scf[0]
    elif n == 2:
        pvtz, pvqz = last_three_scf[0], last_three_scf[1]
    elif n == 3:
        pvdz, pvtz, pvqz = last_three_scf[0], last_three_scf[1], last_three_scf[2]

    # Print a quick diagnostic (optional)
    print(f"[extract_values] {file_path}")
    print(f"  SCF energies kept (oldest->newest): {list(last_three_scf)}")
    print(f"  -> PVDZ: {pvdz}  PVTZ: {pvtz}  PVQZ: {pvqz}  |  Gibbs: {gibbs} | Enthalpy: {enthalpy}")

    return pvdz, pvtz, pvqz, (np.nan if gibbs is None else gibbs), (np.nan if enthalpy is None else enthalpy)


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
        
        pvdz_energy, pvtz_energy, pvqz_energy, gibbs_free_energy, enthalpy = extract_values(file_path)
        
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
                    'Complex Enthalpy Correction': enthalpy,
                })
            elif filename.endswith('Complex_R1.log'):
                data_dict[identification_number].update({
                    'ComplexR1 PVDZ Energy': pvdz_energy,
                    'ComplexR1 PVTZ Energy': pvtz_energy,
                    'ComplexR1 PVQZ Energy': pvqz_energy,
                    'ComplexR1 Gibbs Correction': gibbs_free_energy,
                    'ComplexR1 Enthalpy Correction': enthalpy,
                })
            elif filename.endswith('Complex_R2.log'):
                data_dict[identification_number].update({
                    'ComplexR2 PVDZ Energy': pvdz_energy,
                    'ComplexR2 PVTZ Energy': pvtz_energy,
                    'ComplexR2 PVQZ Energy': pvqz_energy,
                    'ComplexR2 Gibbs Correction': gibbs_free_energy,
                    'ComplexR2 Enthalpy Correction': enthalpy,
                })
            elif filename.endswith('SP.log'):
                data_dict[identification_number].update({
                    'TS PVDZ Energy': pvdz_energy,
                    'TS PVTZ Energy': pvtz_energy,
                    'TS PVQZ Energy': pvqz_energy,
                    'TS Gibbs Correction': gibbs_free_energy,
                    'TS Enthalpy Correction': enthalpy,
                })
            elif filename.endswith('Product.log'):
                data_dict[identification_number].update({
                    'Product PVDZ Energy': pvdz_energy,
                    'Product PVTZ Energy': pvtz_energy,
                    'Product PVQZ Energy': pvqz_energy,
                    'Product Gibbs Correction': gibbs_free_energy,
                    'Product Enthalpy Correction': enthalpy,
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
                    'Complex Enthalpy Correction': enthalpy,
                })
            elif 'IRC' not in filename and filename.endswith('.log'):
                data_dict[identification_number].update({
                    ts_basis_energy_label: pvqz_energy,  # Assign PVQZ energy to TS basis_1 energy
                    'TS Gibbs Correction': gibbs_free_energy,
                    'TS Enthalpy Correction': enthalpy,
                })
            elif filename.endswith('Product.log'):
                data_dict[identification_number].update({
                    product_basis_energy_label: pvqz_energy,  # Assign PVQZ energy to Product basis_1 energy
                    'Product Gibbs Correction': gibbs_free_energy,
                    'Product Enthalpy Correction': enthalpy,
                })
            elif filename.endswith('Complex_R1.log'):
                data_dict[identification_number].update({
                    'ComplexR1 PVQZ Energy': pvqz_energy,  # Assign PVQZ energy to ComplexR1 basis_1 energy
                    'ComplexR1 Gibbs Correction': gibbs_free_energy,
                    'ComplexR1 Enthalpy Correction': enthalpy,
                })
            elif filename.endswith('Complex_R2.log'):
                data_dict[identification_number].update({
                    'ComplexR2 PVQZ Energy': pvqz_energy,  # Assign PVQZ energy to ComplexR2 basis_1 energy
                    'ComplexR2 Gibbs Correction': gibbs_free_energy,
                    'ComplexR2 Enthalpy Correction': enthalpy,
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

df['Enth Separate reagents'] = df['Extrapolated Reagent 1 Energy']+df['Extrapolated Reagent 2 Energy']+df['ComplexR1 Enthalpy Correction']+df['ComplexR2 Enthalpy Correction']
df['Gibbs of separate reagents'] = df['Extrapolated Reagent 1 Energy']+df['Extrapolated Reagent 2 Energy']+df['ComplexR1 Gibbs Correction']+df['ComplexR2 Gibbs Correction']
df['Minimal Enth of separate reagents'] = df['Enth Separate reagents'].min()
df['Minimal Gibbs of separate reagents'] = df['Gibbs of separate reagents'].min()
df['Delta of separate reagents'] = df['Gibbs of separate reagents'] - df['Minimal Gibbs of separate reagents']

df['Gibbs Complex Energy'] = df['Extrapolated Complex Energy'] + df['Complex Gibbs Correction'] 
df['Gibbs TS Energy'] = df['Extrapolated TS Energy'] + df['TS Gibbs Correction'] 
df['Gibbs Product Energy'] = df['Extrapolated Product Energy'] + df['Product Gibbs Correction'] 

df['Enth Complex Energy'] = df['Extrapolated Complex Energy'] + df['Complex Enthalpy Correction']
df['Enth TS Energy'] = df['Extrapolated TS Energy'] + df['TS Enthalpy Correction']
df['Enth Product Energy'] = df['Extrapolated Product Energy'] + df['Product Enthalpy Correction']

df['Separate Reagents'] = 627.5 * (df['Gibbs of separate reagents'] - df['Minimal Gibbs of separate reagents'])
df['Complex Energy'] = 627.5 * (df['Extrapolated Complex Energy'] + df['Complex Gibbs Correction'] - df['Minimal Gibbs of separate reagents'])
df['TS Energy'] = 627.5 * (df['Extrapolated TS Energy'] + df['TS Gibbs Correction'] - df['Minimal Gibbs of separate reagents'])
df['Product Energy'] = 627.5 * (df['Extrapolated Product Energy'] + df['Product Gibbs Correction'] - df['Minimal Gibbs of separate reagents'])

# Ensure 'Complex Energy' is numeric and handle NaN or non-numeric values
df['Separate Reagents'] = pd.to_numeric(df['Separate Reagents'], errors='coerce').fillna(0)
df['Complex Energy'] = pd.to_numeric(df['Complex Energy'], errors='coerce')
df['TS Energy'] = pd.to_numeric(df['TS Energy'], errors='coerce').fillna(0)
df['Product Energy'] = pd.to_numeric(df['Product Energy'], errors='coerce').fillna(0)

df['Enthalpy Separate Reagents'] = 627.5 * (df['Enth Separate reagents'] - df['Enth Separate reagents'].min())
df['Enthalpy Complex Energy'] = 627.5 * (df['Enth Complex Energy'] - df['Enth Separate reagents'].min())
df['Enthalpy TS Energy'] = 627.5 * (df['Enth TS Energy'] - df['Enth Separate reagents'].min())
df['Enthalpy Product Energy'] = 627.5 * (df['Enth Product Energy'] - df['Enth Separate reagents'].min())

# ---------- Kinetics & plotting add-on ----------

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Physical constants (SI)
kB = 1.380649e-23      # J·K^-1
h  = 6.62607015e-34    # J·s
T  = 298.15            # K
Rj = 8.314462618       # J·mol^-1·K^-1
KCAL_TO_J_PER_MOL = 4184.0

# Helper: Eyring prefactor (s^-1)
eyring_prefactor = (kB * T) / h

key_cols = ['Separate Reagents', 'Complex Energy', 'TS Energy', 'Product Energy']
mask_failed = (
    df[key_cols].isna().any(axis=1)
    | (df[['Complex Energy','TS Energy','Product Energy']] == 0).all(axis=1)
)


# --- 1) Reference energies per row for forward direction ---
# If Complex is stable (negative), use Complex as reference; else use Separate Reagents
ref_forward = np.where(df['Complex Energy'].notna() & (df['Complex Energy'] < 0.0),
                       df['Complex Energy'],
                       df['Separate Reagents'])



# Π values for forward selection:
# - For rows with stable complex: use Complex energy in the Boltzmann factor
# - Otherwise: use Separate Reagents
pi_forward_all = pd.Series(
    np.exp(
        - np.where(df['Complex Energy'].notna() & (df['Complex Energy'] < 0.0),
                   df['Complex Energy'],
                   df['Separate Reagents']
        ) * KCAL_TO_J_PER_MOL / (Rj * T)
    ),
    index=df.index
)

# =========================================
# FORWARD: rates + 1dp Product weighting bucket
# =========================================

# --- 2) Special rule: duplicate Separate Reagents (rounded to 2 decimals) ---
# Form groups by SR rounded to 2 decimals (NaNs coerced to 0 which is fine since SR was filled with 0 earlier)
# Form groups by SR rounded to 2 decimals (NaNs coerced to 0 which is fine since SR was filled with 0 earlier)
cr12_sum_5dp = (df['ComplexR1 PVDZ Energy']+df['ComplexR2 PVDZ Energy']).round(5)
df['_CR12_sum_5dp'] = cr12_sum_5dp

# Within each SR group, find the index with the smallest forward barrier TS - SR (NOT TS - Complex),
# because your rule states it’s based on equality of separate reagents.
# ensure mask_failed is the *safe* one and not overwritten
delta_ts_minus_sr_raw = (df['TS Energy'] - df['Separate Reagents']).copy()

# failed rows cannot win
delta_ts_minus_sr_raw[mask_failed] = np.inf
delta_ts_minus_sr_raw[~np.isfinite(delta_ts_minus_sr_raw)] = np.inf

winners = delta_ts_minus_sr_raw.groupby(df['_CR12_sum_5dp']).idxmin()


# Build "effective Π" for percentage: everyone keeps their Π, except
# members of SR-equal groups that are NOT the winner → set Π_eff = 0
pi_eff = pi_forward_all.copy()

for group_val, winner_idx in winners.items():
    members = df.index[(df['_CR12_sum_5dp'] == group_val) & (~mask_failed)]
    if len(members) > 1:
        losers = [i for i in members if i != winner_idx]
        pi_eff.loc[losers] = 0.0

pi_eff_valid = pi_eff.mask(mask_failed, np.nan)
pi_sum = pi_eff_valid.sum(skipna=True)
percentage = pi_eff_valid / pi_sum if pi_sum and np.isfinite(pi_sum) else pd.Series(0.0, index=df.index)


df['Pi (forward ref)'] = pi_forward_all
df['Pi (eff for %)']   = pi_eff
df['Percentage']       = percentage

dg_f_kcal = df['TS Energy'] - ref_forward

# --- 3) Forward rate constants and weighted forward barrier/rate ---
dg_f_Jpermol = dg_f_kcal * KCAL_TO_J_PER_MOL
k_forward = pd.Series(eyring_prefactor * np.exp(-dg_f_Jpermol / (Rj * T)), index=df.index)


# Group-by-1dp rule for SUM of k_forward
df['_CR12_sum_5dp'] = df['Separate Reagents'].round(1)
weighted_k_forward = pd.Series(0.0, index=df.index)
weighted_dgF_kcal  = pd.Series(0.0, index=df.index)

for group_val, group_idxs in df.groupby('_CR12_sum_5dp').groups.items():
    group_idxs = [i for i in group_idxs if not mask_failed.loc[i]]
    if not group_idxs:
        continue
    # Sum forward rates and barrier heights only over valid rows in the 1dp group
    k_sum = k_forward.loc[group_idxs].sum()
    # Winner by (1dp) Percentage among valid rows
    winner_idx = percentage.loc[group_idxs].idxmax()
    weighted_k_forward.loc[winner_idx] = k_sum * percentage.loc[winner_idx]

# Weighted SR
weighted_SR = df['Separate Reagents'] * percentage
weighted_SR_enthalpy = df['Enthalpy Separate Reagents'] * percentage


# =========================================
# REVERSE: rates + 1dp Product weighting bucket
# =========================================

prod_pvdz_5dp = df['Product PVDZ Energy'].round(5)
df['_PROD_5dp'] = prod_pvdz_5dp

delta_ts_minus_prod = (df['TS Energy'] - df['Product Energy']).copy()
delta_ts_minus_prod[mask_failed] = np.inf
delta_ts_minus_prod[~np.isfinite(delta_ts_minus_prod)] = np.inf

winners_rev = delta_ts_minus_prod.groupby(df['_PROD_5dp']).idxmin()

# Base Π for reverse from Product energies
pi_rev = pd.Series(np.exp(- df['Product Energy'] * 4184.0 / (Rj * T)), index=df.index)

# Effective Π for reverse: zero non-winners within duplicate product groups
pi_rev_eff = pi_rev.copy()
for group_val, winner_idx in winners_rev.items():
    members = df.index[(df['_PROD_5dp'] == group_val) & (~mask_failed)]
    if len(members) > 1:
        losers = [i for i in members if i != winner_idx]
        pi_rev_eff.loc[losers] = 0.0

pi_rev_eff_valid = pi_rev_eff.mask(mask_failed, np.nan)
pi_rev_sum = pi_rev_eff_valid.sum(skipna=True)
pct_rev = (pi_rev_eff_valid / pi_rev_sum) if pi_rev_sum and np.isfinite(pi_rev_sum) \
          else pd.Series(0.0, index=df.index)

df['Pi (reverse ref)']   = pi_rev
df['Pi_rev (eff for %)'] = pi_rev_eff
df['Percentage_rev']     = pct_rev


dg_r_kcal = df['TS Energy'] - df['Product Energy']
dg_r_Jpermol = dg_r_kcal * 4184.0
k_reverse = pd.Series(eyring_prefactor * np.exp(-dg_r_Jpermol / (Rj * T)), index=df.index)

df['_PR_1dp'] = df['Product Energy'].round(1)
weighted_k_reverse = pd.Series(0.0, index=df.index)
weighted_dgR_kcal  = pd.Series(0.0, index=df.index)

for group_val, group_idxs in df.groupby('_PR_1dp').groups.items():
    group_idxs = [i for i in group_idxs if not mask_failed.loc[i]]
    if not group_idxs:
        continue
    k_sum = k_reverse.loc[group_idxs].sum()
    winner_idx = pct_rev.loc[group_idxs].idxmax()
    weighted_k_reverse.loc[winner_idx] = k_sum * pct_rev.loc[winner_idx]

# Weighted product
weighted_Product = df['Product Energy'] * pct_rev
weighted_Product_enthalpy = df['Enthalpy Product Energy'] * pct_rev


# --- Store results (failed rows set to NaN so they never contribute to sums/plots) ---
df['Pi (forward ref)']             = pi_forward_all.mask(mask_failed, np.nan)
df['Pi (eff for %)']               = pi_eff.mask(mask_failed, np.nan)
df['Percentage']                   = percentage.mask(mask_failed, np.nan)
df['k_forward (s^-1)']             = k_forward.mask(mask_failed, np.nan)
df['Weighted k_forward (s^-1)']    = weighted_k_forward.mask(mask_failed, np.nan)
df['ΔG‡_forward (kcal/mol)']       = dg_f_kcal.mask(mask_failed, np.nan)
df['Weighted SR']     = weighted_SR.mask(mask_failed, np.nan)
df['Weighted SR (enthalpy)']     = weighted_SR_enthalpy.mask(mask_failed, np.nan)

df['Pi (reverse)']                 = pi_rev.mask(mask_failed, np.nan)
df['Percentage (reverse)']         = pct_rev.mask(mask_failed, np.nan)
df['k_reverse (s^-1)']             = k_reverse.mask(mask_failed, np.nan)
df['Weighted k_reverse (s^-1)']    = (k_reverse * pct_rev).mask(mask_failed, np.nan)
df['ΔG‡_reverse (kcal/mol)']       = dg_r_kcal.mask(mask_failed, np.nan)
df['Weighted Product']             = weighted_Product.mask(mask_failed, np.nan)
df['Weighted Product (enthalpy)']             = weighted_Product_enthalpy.mask(mask_failed, np.nan)

# Optional display columns (percent as 0–100 with two decimals) for Excel
df['Percentage Forward Display'] = (df['Percentage'] * 100).round(2)
df['Percentage Reverse Display'] = (df['Percentage (reverse)'] * 100).round(2)



# ---------- Excel export ----------

# Map to your older naming so your Excel columns are familiar
df['Pi Value Display']     = df['Pi (eff for %)']   # what you show as Pi
df['Percentage Forward']   = df['Percentage']
df['Percentage Reverse']   = df['Percentage (reverse)']
df['Percentage Forward Display'] = (df['Percentage Forward'] * 100).round(2)
df['Percentage Reverse Display'] = (df['Percentage Reverse'] * 100).round(2)
df['Forward Barrier']      = df['ΔG‡_forward (kcal/mol)']
df['Reverse Barrier']      = df['ΔG‡_reverse (kcal/mol)']
df['Forward Rate Constant'] = df['k_forward (s^-1)']
df['Reverse Rate Constant'] = df['k_reverse (s^-1)']


# Optional: clean helper column
df.drop(columns=['_CR12_sum_5dp'], inplace=True)

round_cols = [
    'Separate Reagents','Complex Energy','TS Energy','Product Energy',
    'Pi Value Display','Percentage Forward','Percentage Reverse',
    'Forward Barrier','Reverse Barrier'
]
for col in round_cols:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce').round(2)



preferred = [
    'ID Number',
    'Separate Reagents','Complex Energy','TS Energy','Product Energy',
    'Forward Barrier','Percentage Forward Display','Reverse Barrier','Percentage Reverse Display',
    'Forward Rate Constant','Reverse Rate Constant',
    'Forward Rate (group sum @1dp)', 
    'Weighted k_forward (s^-1)','Weighted k_reverse (s^-1)',
    'ΔG‡_forward (kcal/mol)','ΔG‡_reverse (kcal/mol)',
    'Delta of separate reagents',
    'Enthalpy Separate Reagents','Enthalpy Complex Energy','Enthalpy TS Energy','Enthalpy Product Energy'
    # raw energies / components (include only if present)
    'ComplexR1 PVDZ Energy','ComplexR1 PVTZ Energy','ComplexR1 PVQZ Energy',
    'Extrapolated Reagent 1 Energy','ComplexR1 Gibbs Correction', 'Enth ComplexR1 Enthalpy Correction',
    'ComplexR2 PVDZ Energy','ComplexR2 PVTZ Energy','ComplexR2 PVQZ Energy',
    'Extrapolated Reagent 2 Energy','ComplexR2 Gibbs Correction', 'Enth ComplexR2 Enthalpy Correction',
    'Complex PVDZ Energy','Complex PVTZ Energy','Complex PVQZ Energy',
    'Extrapolated Complex Energy','Complex Gibbs Correction', 'Enth Complex Enthalpy Correction',
    'TS PVDZ Energy','TS PVTZ Energy','TS PVQZ Energy',
    'Extrapolated TS Energy','TS Gibbs Correction', 'Enth TS Enthalpy Correction',
    'Product PVDZ Energy','Product PVTZ Energy','Product PVQZ Energy',
    'Extrapolated Product Energy','Product Gibbs Correction','Enth Product Enthalpy Correction',    
]

ordered = [c for c in preferred if c in df.columns]
#ordered += [c for c in df.columns if c not in ordered]  # append the rest

def eyring_barrier_kcal_from_rate(k):
    # ΔG‡ = - R T ln( (k h) / (kB T) )  [J/mol]  → divide by 4184 for kcal/mol
    if not np.isfinite(k) or k <= 0.0:
        return np.nan
    return -(Rj * T / KCAL_TO_J_PER_MOL) * np.log((k * h) / (kB * T))

# Pull summed weighted rates (safe if columns are missing)
k_f_sum = df['Weighted k_forward (s^-1)'].sum(skipna=True) if 'Weighted k_forward (s^-1)' in df.columns else np.nan
k_r_sum = df['Weighted k_reverse (s^-1)'].sum(skipna=True) if 'Weighted k_reverse (s^-1)' in df.columns else np.nan

# Barriers via Eyring–Polanyi from summed weighted rates
dgF_eyring_kcal = eyring_barrier_kcal_from_rate(k_f_sum)
dgR_eyring_kcal = eyring_barrier_kcal_from_rate(k_r_sum)

# Your requested metric:


# A small Weighted Averages / totals sheet (tweak as you like)
avg_data = pd.DataFrame({
    'Metric': [
        'Sum Weighted k_forward (s^-1)',
        'Weighted ΔG‡_forward from Σ k_forward (kcal/mol)',
        'Sum Weighted k_reverse (s^-1)',
        'Weighted ΔG‡_reverse from Σ k_reverse (kcal/mol)',
        'Mean Reaction Barrier (kcal/mol)',
        'Mean Enthalpy of Reaction (kcal/mol)',

    ],
    'Value': [
        df['Weighted k_forward (s^-1)'].sum(),
        dgF_eyring_kcal,
        df['Weighted k_reverse (s^-1)'].sum(),
        dgR_eyring_kcal,
        df['Weighted Product'].sum() - df['Weighted SR'].sum(),
        df['Weighted Product (enthalpy)'].sum() - df['Weighted SR (enthalpy)'].sum(),
    ]
})

# Columns to label as "Calc failed" in the spreadsheet view
cols_to_flag = [
    'Separate Reagents','Complex Energy','TS Energy','Product Energy',
    'Forward Barrier','Reverse Barrier',
    'Percentage Forward Display','Percentage Reverse Display',
    'Forward Rate Constant','Reverse Rate Constant',
    'Weighted forward barrier (per 1dp rule)',
    'Weighted k_forward (s^-1)','Weighted k_reverse (s^-1)',
    'ΔG‡_forward (kcal/mol)','ΔG‡_reverse (kcal/mol)',
    'Weighted ΔG‡_forward (kcal)','Weighted ΔG‡_reverse (kcal)'
]
# map to your actual column names if you use different ones:
name_map = {
    'Forward Barrier': 'ΔG‡_forward (kcal/mol)',
    'Reverse Barrier': 'ΔG‡_reverse (kcal/mol)',
    'Forward Rate Constant': 'k_forward (s^-1)',
    'Reverse Rate Constant': 'k_reverse (s^-1)'
}
cols_to_flag = [name_map.get(c, c) for c in cols_to_flag]
cols_to_flag = [c for c in cols_to_flag if c in df.columns]

df = df.sort_values(by='ID Number', ascending=True, na_position='last')

df_excel = df.copy()
df_excel.loc[mask_failed, cols_to_flag] = "Calc failed"

with pd.ExcelWriter('BOLTCAR_results.xlsx') as writer:
    # Ensure 'ID Number' is a normal column (not index)
    df_reset = df.reset_index(drop=True)
    df_reset.to_excel(writer, index=False, sheet_name='Full Results', columns=ordered)
    avg_data.to_excel(writer, index=False, sheet_name='Weighted Averages')

print("Data extraction complete. The results are saved in 'BOLTCAR_results.xlsx'.")



# ================== PLOTTING ==================


# Filter to relevant rows (Pi > 0 / Percentage > 0); tweak if you prefer reverse as well
df_filtered = df[df['Percentage'] > 0].copy()
if 'ID Number' not in df_filtered.columns:
    df_filtered['ID Number'] = df_filtered.index.astype(str)
df_filtered['ID Number'] = df_filtered['ID Number'].astype(str)

# Percent displays (if not already created)
if 'Percentage Forward Display' not in df_filtered.columns:
    df_filtered['Percentage Forward Display'] = (df_filtered.get('Percentage Forward', df_filtered['Percentage']) * 100).round(2)
if 'Percentage Reverse Display' not in df_filtered.columns:
    df_filtered['Percentage Reverse Display'] = (df_filtered.get('Percentage Reverse', df_filtered.get('Percentage (reverse)', 0.0)) * 100).round(2)

# Energies for y-limits (ignore NaNs)
energies_stack = np.vstack([
    df_filtered['Separate Reagents'].to_numpy(dtype=float),
    df_filtered['Complex Energy'].to_numpy(dtype=float, copy=True),
    df_filtered['TS Energy'].to_numpy(dtype=float),
    df_filtered['Product Energy'].to_numpy(dtype=float)
])
finite_vals = energies_stack[np.isfinite(energies_stack)]
if finite_vals.size == 0:
    min_energy, max_energy = -1.0, 1.0
else:
    pad = max(1.0, 0.05 * (finite_vals.max() - finite_vals.min()))
    min_energy = finite_vals.min() - pad
    max_energy = finite_vals.max() + pad
df_valid = df[~mask_failed].copy()

# Example: totals for legends
sum_weighted_kf = float(df_valid.get('Weighted k_forward (s^-1)', pd.Series(dtype=float)).sum())
sum_weighted_kr = float(df_valid.get('Weighted k_reverse (s^-1)', pd.Series(dtype=float)).sum())
label_forward = f"Σ Weighted kf = {sum_weighted_kf:.2e}" if np.isfinite(sum_weighted_kf) else ""
label_reverse = f"Σ Weighted kr = {sum_weighted_kr:.2e}" if np.isfinite(sum_weighted_kr) else ""

# Your plotting filter should use df_valid (and e.g., df_filtered = df_valid[df_valid['Percentage'] > 0])
df_filtered = df_valid[df_valid['Percentage'] > 0].copy()

# Output directory
output_dir = os.path.abspath("plots")
os.makedirs(output_dir, exist_ok=True)


# ---------- Overall scatter plots (all IDs) ----------

# Separate Reagents (optional but handy)
plt.figure(figsize=(10, 5))
plt.scatter(df_filtered['ID Number'], df_filtered['Separate Reagents'], label='Separate Reagents')
plt.xlabel('ID Number')
plt.ylabel('Energy (kcal/mol)')
plt.title('Separate Reagents Energies (Pi > 0)')
plt.xticks(rotation=90)
plt.legend()
plt.ylim(min_energy, max_energy)
plt.grid(True, which='both', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "sr_energies.png"), dpi=300, bbox_inches='tight')
plt.close()

# Complex
plt.figure(figsize=(10, 5))
plt.scatter(df_filtered['ID Number'], df_filtered['Complex Energy'],
            color='b',
            label=f'Complex Energy\n{label_forward}' if label_forward else 'Complex Energy')
plt.xlabel('ID Number')
plt.ylabel('Energy (kcal/mol)')
plt.title('Complex Energies (Pi > 0)')
plt.xticks(rotation=90)
plt.legend()
plt.ylim(min_energy, max_energy)
plt.grid(True, which='both', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "complex_energies.png"), dpi=300, bbox_inches='tight')
plt.close()

# TS
plt.figure(figsize=(10, 5))
plt.scatter(df_filtered['ID Number'], df_filtered['TS Energy'], label='TS Energy')
plt.xlabel('ID Number')
plt.ylabel('Energy (kcal/mol)')
plt.title('Transition State Energies (Pi > 0)')
plt.xticks(rotation=90)
plt.legend()
plt.ylim(min_energy, max_energy)
plt.grid(True, which='both', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "ts_energies.png"), dpi=300, bbox_inches='tight')
plt.close()

# Product
plt.figure(figsize=(10, 5))
plt.scatter(df_filtered['ID Number'], df_filtered['Product Energy'],
            color='g',
            label=f'Product Energy\n{label_reverse}' if label_reverse else 'Product Energy')
plt.xlabel('ID Number')
plt.ylabel('Energy (kcal/mol)')
plt.title('Product Energies (Pi > 0)')
plt.xticks(rotation=90)
plt.legend()
plt.ylim(min_energy, max_energy)
plt.grid(True, which='both', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "product_energies.png"), dpi=300, bbox_inches='tight')
plt.close()

# Percentages (use the display %)
plt.figure(figsize=(10, 5))
plt.scatter(df_filtered['ID Number'], df_filtered['Percentage Forward Display'], label='Forward %')
if 'Percentage Reverse Display' in df_filtered:
    plt.scatter(df_filtered['ID Number'], df_filtered['Percentage Reverse Display'], label='Reverse %')
plt.xlabel('ID Number')
plt.ylabel('Percentage (%)')
plt.title('Percentages (Pi > 0)')
plt.xticks(rotation=90)
plt.legend()
plt.ylim(-1, 101)
plt.grid(True, which='both', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "percentages.png"), dpi=300, bbox_inches='tight')
plt.close()

# ---------- 2x2 combined page (whole dataset) ----------
fig, axs = plt.subplots(2, 2, figsize=(16, 12))

axs[0, 0].scatter(df_filtered['ID Number'], df_filtered['Complex Energy'], color="#24BB7A", label=f'Complex Energy\n{label_forward}')
axs[0, 0].set_title('Complex Energies')
axs[0, 0].set_xlabel('ID Number')
axs[0, 0].set_ylabel('Energy (kcal/mol)')
axs[0, 0].set_ylim(min_energy, max_energy)
axs[0, 0].legend(); axs[0, 0].grid(True, linestyle='--', alpha=0.4); axs[0, 0].tick_params(axis='x', rotation=90)

axs[0, 1].scatter(df_filtered['ID Number'], df_filtered['TS Energy'], color="#065143", label='TS Energy')
axs[0, 1].set_title('Transition State Energies')
axs[0, 1].set_xlabel('ID Number')
axs[0, 1].set_ylabel('Energy (kcal/mol)')
axs[0, 1].set_ylim(min_energy, max_energy)
axs[0, 1].legend(); axs[0, 1].grid(True, linestyle='--', alpha=0.4); axs[0, 1].tick_params(axis='x', rotation=90)

axs[1, 0].scatter(df_filtered['ID Number'], df_filtered['Product Energy'], color="#F38503", label=f'Product Energy\n{label_reverse}')
axs[1, 0].set_title('Product Energies')
axs[1, 0].set_xlabel('ID Number')
axs[1, 0].set_ylabel('Energy (kcal/mol)')
axs[1, 0].set_ylim(min_energy, max_energy)
axs[1, 0].legend(); axs[1, 0].grid(True, linestyle='--', alpha=0.4); axs[1, 0].tick_params(axis='x', rotation=90)

axs[1, 1].scatter(df_filtered['ID Number'], df_filtered['Percentage Forward Display'], color="#24BB7A", label='Forward %')
if 'Percentage Reverse Display' in df_filtered:
    axs[1, 1].scatter(df_filtered['ID Number'], df_filtered['Percentage Reverse Display'], color="#F38503", label='Reverse %')
axs[1, 1].set_title('Reaction Percentages')
axs[1, 1].set_xlabel('ID Number')
axs[1, 1].set_ylabel('Percentage (%)')
axs[1, 1].set_ylim(-1, 101)
axs[1, 1].legend(); axs[1, 1].grid(True, linestyle='--', alpha=0.4); axs[1, 1].tick_params(axis='x', rotation=90)

plt.tight_layout()
pdf_path_combined = os.path.join(output_dir, "BOLTCAR_combined_plots.pdf")
plt.savefig(pdf_path_combined, dpi=300, bbox_inches='tight')
plt.close()
print(f"Combined 2x2 plots saved as '{pdf_path_combined}'.")

# ---------- Per-ID bar plots into a single PDF ----------
per_id_pdf = os.path.join(output_dir, "BOLTCAR_perID_plots.pdf")
with PdfPages(per_id_pdf) as pdf:
    for _, row in df_filtered.sort_values('ID Number').iterrows():
        fig, ax = plt.subplots(figsize=(8.5, 5.5))
        labels = ['Separate Reagents', 'Complex', 'TS', 'Product']
        vals = [
            float(row.get('Separate Reagents', np.nan)),
            float(row.get('Complex Energy', np.nan)),
            float(row.get('TS Energy', np.nan)),
            float(row.get('Product Energy', np.nan))
        ]
        ax.bar(labels, vals)
        pct_f = float(row.get('Percentage Forward Display', np.nan))
        pct_r = float(row.get('Percentage Reverse Display', np.nan))
        ttl = f"ID {row['ID Number']} | %Fwd={pct_f:.2f}%"
        if np.isfinite(pct_r):
            ttl += f" | %Rev={pct_r:.2f}%"
        ax.set_title(ttl)
        ax.set_ylabel('Energy (kcal/mol)')
        ax.set_ylim(min_energy, max_energy)
        ax.grid(axis='y', linestyle='--', alpha=0.4)
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

print(f"Per-ID bar plots written to '{per_id_pdf}'.")
print(f"All plots saved in '{output_dir}'.")


