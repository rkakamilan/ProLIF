# ProLIF Codebase Architecture

## Directory Structure
```
prolif/
├── __init__.py           # Main package exports
├── _version.py          # Version information
├── fingerprint.py       # Core Fingerprint class
├── molecule.py          # Molecule wrapper class
├── residue.py           # Residue identification system
├── ifp.py              # Interaction fingerprint data structure
├── utils.py            # Utility functions
├── exceptions.py       # Custom exceptions
├── interactions/       # Interaction detection algorithms
├── plotting/           # Visualization modules
└── data/              # Example data files
```

## Core Components

### 1. Fingerprint Class (`prolif.fingerprint.Fingerprint`)
- **Purpose**: Main class for calculating interaction fingerprints
- **Features**: 
  - Runs analysis on single structures or docking results
  - Supports parallel processing with multiprocessing
  - Stores results as IFP (Interaction Fingerprint) dictionaries
- **Key Methods**:
  - `generate()`: Generate fingerprint between two molecules
  - `run_from_iterable()`: Process multiple ligands against one protein
  - `to_dataframe()`, `to_bitvectors()`: Convert results to different formats

### 2. Interactions (`prolif.interactions`)
- **Base classes**: `interactions/base.py` (Distance, SingleAngle, DoubleAngle)
- **Implementations**: `interactions/interactions.py` (HBDonor, HBAcceptor, PiStacking, etc.)
- **Special interactions**: `interactions/water_bridge.py` (bridged interactions)

### 3. Molecule (`prolif.molecule.Molecule`)
- **Purpose**: Wrapper around RDKit molecules
- **Key Features**:
  - Handles conversion from various formats
  - Provides geometric properties (centroid, xyz coordinates)
  - **Note**: No longer supports MDAnalysis integration

### 4. Residue System (`prolif.residue.ResidueId`)
- **Purpose**: Unique identifiers for protein/ligand residues
- **Features**: Handles chain, residue number, and insertion code parsing

### 5. IFP Data Structure (`prolif.ifp.IFP`)
- **Purpose**: Storage for interaction fingerprints
- **Features**:
  - Dictionary-like interface mapping residue pairs to interaction data
  - Supports conversion to DataFrames and bitvectors

## Data Flow
1. **Input**: PDB/SDF files → Molecule objects (via RDKit)
2. **Processing**: Molecule objects → Interaction detection (via interaction classes)
3. **Storage**: Detected interactions → IFP storage
4. **Output**: IFP → Analysis outputs (DataFrames, bitvectors, plots)

## Key Design Patterns

### MDAnalysis Removal Strategy
- **Conditional imports**: Graceful handling of missing MDAnalysis
- **Dummy classes**: Placeholder classes to maintain API compatibility
- **Error messages**: Informative errors with migration guidance
- **Feature deprecation**: Clear marking of unsupported functionality

### Interaction Detection
- **Plugin architecture**: Modular interaction classes
- **Base classes**: Common patterns (Distance, SingleAngle, DoubleAngle)
- **Metadata support**: Rich interaction information beyond simple detection

### Parallel Processing
- **Multiprocessing support**: Built-in parallel execution for large datasets
- **Progress tracking**: tqdm integration for user feedback
- **Memory management**: Efficient handling of large molecular datasets

## Testing Structure
- **Main tests**: `tests/` directory
- **Plotting tests**: `tests/plotting/`
- **Fixtures**: `conftest.py` with example molecules
- **Data**: Test fixtures from `prolif/data/`

## Important Notes
- **RDKit-centric**: Primary dependency is RDKit for chemical informatics
- **Docking focus**: Optimized for docking results rather than MD trajectories
- **Licensing**: Apache 2.0 license enables broader usage
- **Python 3.10+**: Modern Python version requirements