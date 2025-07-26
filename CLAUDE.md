# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 【MUST GLOBAL】Gemini活用（プロジェクトのCLAUDE.mdより優先）

### 三位一体の開発原則
ユーザーの**意思決定**、Claudeの**分析と実行**、Geminiの**検証と助言**を組み合わせ、開発の質と速度を最大化する：
- **ユーザー**：プロジェクトの目的・要件・最終ゴールを定義し、最終的な意思決定を行う**意思決定者**
  - 反面、具体的なコーディングや詳細な計画を立てる力、タスク管理能力ははありません。
- **Claude**：高度な計画力・高品質な実装・リファクタリング・ファイル操作・タスク管理を担う**実行者**
  - 指示に対して忠実に、順序立てて実行する能力はありますが、意志がなく、思い込みは勘違いも多く、思考力は少し劣ります。
- **Gemini**：深いコード理解・Web検索 (Google検索) による最新情報へのアクセス・多角的な視点からの助言・技術的検証を行う**助言者**
  - プロジェクトのコードと、インターネット上の膨大な情報を整理し、的確な助言を与えてくれますが、実行力はありません。

### 実践ガイド
- **ユーザーの要求を受けたら即座に`gemini -p <質問内容>`で壁打ち**を必ず実施
- Geminiの意見を鵜呑みにせず、1意見として判断。聞き方を変えて多角的な意見を抽出
- Claude Code内蔵のWebSearchツールは使用しない
- Geminiがエラーの場合は、聞き方を工夫してリトライ：
  - ファイル名や実行コマンドを渡す（Geminiがコマンドを実行可能）
  - 複数回に分割して聞く

### 主要な活用場面
1. **実現不可能な依頼**: Claude Codeでは実現できない要求への対処 (例: `今日の天気は？`)
2. **前提確認**: ユーザー、Claude自身に思い込みや勘違い、過信がないかどうか逐一確認 (例: `この前提は正しいか？`）
3. **技術調査**: 最新情報・エラー解決・ドキュメント検索・調査方法の確認（例: `Rails 7.2の新機能を調べて`）
4. **設計検証**: アーキテクチャ・実装方針の妥当性確認（例: `この設計パターンは適切か？`）
5. **問題解決**: Claude自身が自力でエラーを解決できない場合に対処方法を確認 (例: `この問題の対処方法は？`)
6. **コードレビュー**: 品質・保守性・パフォーマンスの評価（例: `このコードの改善点は？`）
7. **計画立案**: タスクの実行計画レビュー・改善提案（例: `この実装計画の問題点は？`）
8. **技術選定**: ライブラリ・手法の比較検討 （例: `このライブラリは他と比べてどうか？`

* [Ref](https://zenn.dev/tksfjt1024/articles/5e88385bfb69fd)


## Project Overview

ProLIF (Protein-Ligand Interaction Fingerprints) is a Python package for generating interaction fingerprints between proteins, ligands, DNA, and RNA molecules from docking simulations and experimental structures. 

**⚠️ IMPORTANT: MDAnalysis Dependency Removed**
As of the latest version, MDAnalysis dependency has been removed to resolve LGPL3 licensing issues. The package now focuses on docking results analysis (PDB/PDBQT proteins + SDF ligands) and no longer supports MD trajectory analysis. This change enables use under more permissive licenses.

## Commands

### Development and Testing
- `poe format` - Format code using ruff and black
- `poe lint` - Lint code using ruff 
- `poe style-fix` - Fix formatting and linting issues
- `poe type-check` - Run mypy type checking
- `poe test` - Run test suite with pytest
- `poe docs` - Build documentation with Sphinx
- `poe check` - Run all checks (style, type, test, docs)

### Style Checks Only
- `poe format-check` - Check if code needs formatting
- `poe lint-check` - Check if code needs linting
- `poe style-check` - Run both formatting and linting checks

### Single Test
```bash
pytest tests/test_fingerprint.py::TestFingerprint::test_method_name
```

### Dependencies
The project uses uv for dependency management. Install with:
```bash
uv sync
```

## Architecture

### Core Components

1. **Fingerprint (`prolif.fingerprint.Fingerprint`)** - Main class for calculating interaction fingerprints
   - Runs analysis on MD trajectories or single structures
   - Supports parallel processing with multiprocessing
   - Stores results as IFP (Interaction Fingerprint) dictionaries

2. **Interactions (`prolif.interactions`)** - Interaction detection algorithms
   - Base classes in `interactions/base.py` define common patterns (Distance, SingleAngle, DoubleAngle)
   - Concrete implementations in `interactions/interactions.py` (HBDonor, HBAcceptor, PiStacking, etc.)
   - Special water bridge interactions in `interactions/water_bridge.py`

3. **Molecule (`prolif.molecule.Molecule`)** - Wrapper around RDKit molecules with MDAnalysis integration
   - Handles conversion between MDAnalysis AtomGroups and RDKit Mol objects
   - Provides geometric properties (centroid, xyz coordinates)

4. **Residue (`prolif.residue.ResidueId`)** - Residue identification system
   - Unique identifiers for protein/ligand residues
   - Handles chain, residue number, and insertion code parsing

5. **IFP (`prolif.ifp.IFP`)** - Data structure for storing interaction fingerprints
   - Dictionary-like interface mapping residue pairs to interaction data
   - Supports conversion to DataFrames and bitvectors

### Data Flow

1. Input structures (PDB/SDF files) → Molecule objects (via RDKit)
2. Molecule objects → Interaction detection (via interaction classes)
3. Detected interactions → IFP storage
4. IFP → Analysis outputs (DataFrames, bitvectors, plots)

### Supported Input Formats

✅ **Supported** (Docking Analysis):
- PDB files: Protein structures
- SDF files: Ligand structures from docking
- MOL2 files: Limited support via RDKit

❌ **No Longer Supported** (MD Trajectory Analysis):
- PDBQT files: Requires MDAnalysis conversion
- Trajectory formats (XTC, TRR, DCD): MD simulation analysis
- `Molecule.from_mda()`: Use `Molecule.from_rdkit()` instead
- `select_over_trajectory()`: Use RDKit molecule methods

### Key Dependencies

- **RDKit**: Chemical informatics and molecular representations (primary)
- **NetworkX**: Network analysis for interaction networks
- **Pandas**: Data manipulation and analysis
- **NumPy/SciPy**: Numerical computations

### Removed Dependencies

- **MDAnalysis**: ❌ Removed due to LGPL3 licensing. Functions requiring MDAnalysis now raise `NotImplementedError` with migration guidance.

## Migration Guide (MDAnalysis Removal)

### Code Changes Required

**Old (MDAnalysis-based):**
```python
import MDAnalysis as mda
import prolif

# Load from trajectory
u = mda.Universe("protein.pdb", "trajectory.xtc")
ligand_ag = u.select_atoms("resname LIG")
protein_ag = u.select_atoms("protein")

ligand = prolif.Molecule.from_mda(ligand_ag)
protein = prolif.Molecule.from_mda(protein_ag)
```

**New (RDKit-based):**
```python
from rdkit import Chem
import prolif

# Load from static structures
protein_mol = Chem.MolFromPDBFile("protein.pdb", removeHs=False)
protein = prolif.Molecule.from_rdkit(protein_mol)

# Load ligands from SDF
ligands = prolif.sdf_supplier("ligands.sdf")
```

### Deprecated Functions

| Function | Status | Alternative |
|----------|--------|------------|
| `Molecule.from_mda()` | ❌ Removed | Use `Molecule.from_rdkit()` |
| `pdbqt_supplier` | ❌ Removed | Convert PDBQT to PDB/SDF first |
| `select_over_trajectory()` | ❌ Removed | Use RDKit molecule methods |
| Water bridge interactions | ⚠️ Simplified | Basic implementation without MD features |

### Error Messages

All removed functions provide clear error messages with migration guidance.

### Testing Structure

- Main tests in `tests/` directory
- Plotting tests in `tests/plotting/`
- Test fixtures defined in `conftest.py` with example molecules
- Uses pytest with molecular data fixtures from `prolif/data/`

### Data Files

- Example molecular data in `prolif/data/` (MOL2, PDB, XTC files)
- Vina docking outputs in `prolif/data/vina/`
- Network visualization assets in `prolif/plotting/network/`

## Code Style

- Uses ruff for linting and formatting
- Uses black for Jupyter notebook formatting  
- Type hints enforced with mypy
- Follows PEP 8 style guidelines
- No use of `any` or `unknown` types in TypeScript (though this is a Python project)