# BOLTCAR

**Boltzmann Optimized Likelihood Transition-state Calculation and Analysis for Reactions**

BOLTCAR is a high-throughput Python-based pipeline to automatically identify, refine, and analyze transition states (TS) and reaction pathways based on Gaussian output files. It extends the FASTCAR approach by integrating conformer pruning (CREST + RMSD), intrinsic reaction coordinate (IRC) validation, stationary point analysis, and Boltzmann-weighted kinetic modeling.

---

## 🛠 Installation and Setup

Clone the repository:

```bash
git clone https://github.com/Lise345/BOLTCAR.git
cd BOLTCAR
```

Ensure you have:
- Python 3.7+
- Gaussian 16
- CREST, xtb, sPyRMSD, OpenBabel
- SLURM scheduler (for job submission)
- Required Python packages (`pandas`, `numpy`, `matplotlib`, `openpyxl`)

> ⚠️ **Important**: The `.sub` submission scripts include SLURM `module load` commands. **You must adapt these modules to match your cluster environment**, e.g., `module load Gaussian/G16`, `module load xtb`, etc.

---

## 🚀 How to Use

### 1. Prepare Input

You need a Gaussian `.out` file for a **transition state** geometry (with at least one imaginary frequency):

```bash
python Script_1_CREST_INPUT.py your_TS_file.out
```

This will:
- Extract geometry and frequency
- Update `parameters.txt` with inferred properties
- Create and launch `CREST` and pruning jobs via `Script_1_RMSD_INPUT.py`

Each script then triggers the next via `.sub` files and job dependencies.

---

### 2. Workflow Overview

| Step | Script                          | Description                                                          |
|------|---------------------------------|----------------------------------------------------------------------|
| 1    | `Script_1_CREST_INPUT.py`       | Extracts info, runs CREST, submits conformer search jobs            |
| 2    | `Script_1_RMSD_INPUT.py`        | Prunes conformers via RMSD, submits TS optimizations                |
| 3    | `Script_2_IRC_Calculation.py`   | Validates TS via IRC, filters on frequency, energy, and redundancy |
| 4    | `Script_3_StationaryPoints.py`  | Distinguishes product and complex, submits SP calculations          |
| 5    | `Script_4_Reactants.py`         | Extracts separate reagents and launches optimization jobs           |
| 6    | `Script_5_Results_BOLTCAR.py`   | Collects results, calculates rate constants, generates plots        |

---

## 🧾 Input File

You must provide a Gaussian `.out` file with:
- A complete frequency calculation
- Geometry and standard orientation
- Charge and multiplicity

---

## ⚙️ Configuration: `parameters.txt`

This file is used throughout the pipeline to control all major inputs and SLURM parameters.

---

### CREST Parameters

| Parameter        | Example     | Description                                      |
|------------------|-------------|--------------------------------------------------|
| `CREST version`  | `default`   | Select CREST version (`default`, `3.0`, etc.)    |
| `CREST solvent`  | `none`      | Add if CREST solvation is needed (e.g. `--alpb`) |

---

### RMSD Pruning Parameters

| Parameter          | Example | Description                                                       |
|--------------------|---------|-------------------------------------------------------------------|
| `RMSD threshold`   | `0.5`   | Maximum RMSD (Å) to consider TS geometries as equivalent          |
| `Energy threshold` | `0.05`  | Energy diff in kcal/mol below which two TSs are considered similar |
| `Energy window`    | `10`    | Maximum energy range considered in pruning (kcal/mol)            |
| `B_threshold`      | `2`     | Rotational constant RMSD threshold in % for additional filtering |

---

### DFT Parameters

| Parameter                  | Example                    | Description                                                                 |
|----------------------------|----------------------------|-----------------------------------------------------------------------------|
| `Functional`               | `m062x`                    | Functional to use                                                           |
| `Dispersion`              | `empiricaldispersion=gd3` | Dispersion method (or `none`)                                               |
| `Basis`                    | `cbs`                      | Basis set (`cbs` triggers cc-pVDZ/VTZ/VQZ extrapolation)                    |
| `DFT solvent`              | `none`                     | Solvent model (e.g., `scrf=(smd,solvent=chloroform)`) or `none`             |
| `Charge`                   | `0`                        | Total molecular charge                                                      |
| `Multiplicity`             | `1`                        | Spin multiplicity                                                           |
| `Time for TS calcs`        | `40`                       | SLURM job time in hours for TS optimization                                 |
| `Time for IRC calcs`       | `40`                       | SLURM job time in hours for IRC calcs                                       |
| `Time for stationary calcs`| `40`                       | SLURM job time in hours for stationary point calcs                          |

---

### Geometry Mapping

| Parameter           | Example                              | Description                                                           |
|---------------------|--------------------------------------|-----------------------------------------------------------------------|
| `molecule1_atoms`   | `3 8 9 12 13 14 17 ...`              | Atom indices (space-separated) for molecule 1 (reagent 1)            |
| `molecule2_atoms`   | `1 2 4 5 6 7 ...`                    | Atom indices for molecule 2 (reagent 2)                              |
| `atom1`             | `1`                                  | Atom index in molecule 1 forming the bond                            |
| `atom2`             | `12`                                 | Atom index in molecule 2 forming the bond                            |
| `size_molecule`     | `45`                                 | Total number of atoms in the molecule                                |
| `CC1_in`            | `"1 12"`                             | The atoms that define the C–C bond used for IRC direction            |

---

### Used Paths

| Parameter   | Example                                                                 |
|------------|-------------------------------------------------------------------------|
| `rootdir`  | `/rhea/scratch/.../BOLTCAR`                                             |
| `bin`      | `/vscmnt/.../bin/`                                                      |

---

## 📂 Output

- `BOLTCAR_results.xlsx`: rate constants, Boltzmann weights, energies
- `plots/`: comparison plots (PNG + combined PDF)
- `TS_analysis.xlsx`: intermediate summary tables

---

## 📌 Notes

- You must run on a SLURM-based HPC system.
- Only the first script is launched manually. The rest trigger through SLURM dependencies.
- All `.sub` scripts contain `module load` statements — change these for your cluster.

---

## 📧 Contact

If you use this tool in a publication, please cite it appropriately. Feel free to open an issue or contact the repository maintainer for questions or feedback.
