# ProLIF Project Overview

## Project Purpose
ProLIF (Protein-Ligand Interaction Fingerprints) is a Python package for generating interaction fingerprints between proteins, ligands, DNA, and RNA molecules from docking simulations and experimental structures.

**Key Focus**: The package has been modified to be an "MDAnalysis-free version" that focuses on docking results analysis only, with the MDAnalysis dependency removed to resolve LGPL3 licensing issues.

## Current Package Name
- **Package name**: `prolif-lite` (renamed from original `prolif`)
- **Description**: Interaction Fingerprints for protein-ligand docking analysis (MDAnalysis-free version)
- **Version**: 2.1.3 (as of current codebase)
- **License**: Apache License 2.0

## Supported Features
✅ **Currently Supported**:
- Protein-ligand interaction analysis from docking results
- PDB protein structure support  
- SDF ligand structure support
- MOL2 file support (limited via RDKit)
- All core interaction types (hydrogen bonds, hydrophobic, π-stacking, etc.)
- Interaction fingerprint generation and analysis
- Export to DataFrames and bitvectors

❌ **No Longer Supported** (due to MDAnalysis removal):
- MD trajectory analysis (XTC, TRR, DCD formats)
- PDBQT file support (requires conversion to PDB/SDF)
- `Molecule.from_mda()` method
- `select_over_trajectory()` function

## Key Dependencies
- **RDKit**: Chemical informatics and molecular representations (primary)
- **NetworkX**: Network analysis for interaction networks
- **Pandas**: Data manipulation and analysis
- **NumPy/SciPy**: Numerical computations
- **Removed**: MDAnalysis (due to LGPL3 licensing issues)

## Development Focus
This is a specialized version focused on docking analysis rather than MD trajectory analysis, enabling use under more permissive licenses.