# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 【MUST GLOBAL】Gemini活用（プロジェクトのCLAUDE.mdより優先）

### 三位一体の開発原則
ユーザーの**意思決定**、Claudeの**分析と実行**、Geminiの**検証と助言**を組み合わせ、開発の質と速度を最大化する：
- **ユーザー**：プロジェクトの目的・要件・最終ゴールを定義し、最終的な意思決定を行う**意思決定者**
  - 反面、具体的なコーディングや詳細な計画を立てる力、タスク管理能力はありません。
- **Claude**：高度な計画力・高品質な実装・リファクタリング・ファイル操作・タスク管理を担う**実行者**
  - 指示に対して忠実に、順序立てて実行する能力はありますが、意志がなく、思い込みや勘違いも多く、思考力は少し劣ります。
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
2. **前提確認**: ユーザー、Claude自身に思い込みや勘違い、過信がないかどうか逐一確認 (例: `この前提は正しいか？`)
3. **技術調査**: 最新情報・エラー解決・ドキュメント検索・調査方法の確認（例: `Rails 7.2の新機能を調べて`）
4. **設計検証**: アーキテクチャ・実装方針の妥当性確認（例: `この設計パターンは適切か？`）
5. **問題解決**: Claude自身が自力でエラーを解決できない場合に対処方法を確認 (例: `この問題の対処方法は？`)
6. **コードレビュー**: 品質・保守性・パフォーマンスの評価（例: `このコードの改善点は？`）
7. **計画立案**: タスクの実行計画レビュー・改善提案（例: `この実装計画の問題点は？`）
8. **技術選定**: ライブラリ・手法の比較検討（例: `このライブラリは他と比べてどうか？`）

* [Ref](https://zenn.dev/tksfjt1024/articles/5e88385bfb69fd)


## Project Overview

**prolif-lite** (v2.1.3) は ProLIF (Protein-Ligand Interaction Fingerprints) の MDAnalysis 除去フォーク。  
ドッキングシミュレーション結果（PDB/SDF）からタンパク質-リガンド相互作用フィンガープリントを生成する Python パッケージ。

**フォーク元**: https://github.com/chemosim-lab/ProLIF  
**フォークの目的**: LGPL3 ライセンスを持つ MDAnalysis を除去し、Apache License 2.0 のみで配布可能にする。

**⚠️ IMPORTANT: MDAnalysis Dependency Removed**  
MDAnalysis 依存を除去済み。MD軌跡解析は非対応。ドッキング結果解析（PDB タンパク質 + SDF リガンド）にフォーカス。

## Release & PyPI Publishing

### リリース手順
1. `prolif/_version.py` のバージョンを更新してコミット・push
2. GitHub で Release を作成（タグ例: `v2.1.4`）
3. `build.yml` が自動起動し、ビルド後に PyPI へ自動 publish される

### PyPI Trusted Publishing
- `PYPI_TOKEN` シークレットは**不要**（OIDC Trusted Publishing を使用）
- 新しいリポジトリで初めて publish する場合のみ、PyPI 側で一度だけ設定が必要:
  1. https://pypi.org/manage/project/prolif-lite/ を開く
  2. 左メニュー **「Publishing」** タブ → **「Add a new publisher」**
  3. GitHub Actions を選択し以下を入力:
     - Owner: `rkakamilan` / Repository: `ProLIF` / Workflow: `build.yml`

### バージョニング方針
- upstream (chemosim-lab/ProLIF) との連続性を重視し patch bump を基本とする
- 新しい公開 API の追加や動作変更を伴う場合は minor bump

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

1. **Fingerprint (`prolif.fingerprint.Fingerprint`)** - 相互作用フィンガープリントの計算クラス
   - `run_from_iterable()` で複数のリガンドポーズを解析（ドッキング用途）
   - `run()` は MDAnalysis 依存のため NotImplementedError を送出（スタブとして残存）
   - 並列処理対応（multiprocessing）。`PROLIF_N_JOBS` 環境変数で並列数を制御可能
   - 結果を IFP (Interaction Fingerprint) 辞書として保持

2. **Interactions (`prolif.interactions`)** - 相互作用検出アルゴリズム
   - 基底クラス: `interactions/base.py`（Distance, SingleAngle, DoubleAngle, BasePiStacking）
   - 具体実装: `interactions/interactions.py`（HBDonor, HBAcceptor, PiStacking 等）
   - 水分子仲介相互作用: `interactions/water_bridge.py`（`run_from_iterable()` のみ対応）

3. **Molecule (`prolif.molecule.Molecule`)** - RDKit Mol のラッパー
   - `from_rdkit()` で RDKit Mol から生成（主要インタフェース）
   - `from_mda()` は NotImplementedError（MDAnalysis 除去済み）
   - 幾何学プロパティ（centroid, xyz 座標）を提供

4. **Residue (`prolif.residue.ResidueId`)** - 残基識別システム
   - タンパク質/リガンド残基のユニーク識別子
   - チェーン・残基番号・挿入コードのパース

5. **IFP (`prolif.ifp.IFP`)** - 相互作用フィンガープリント格納データ構造
   - 残基ペアを相互作用データにマッピングする辞書ライクインタフェース
   - DataFrame・bitvector への変換対応
   - `IFP.interactions()` イテレータで全相互作用データをフラットに取得可能

### Data Flow

1. PDB/SDF ファイル → Molecule オブジェクト（via RDKit）
2. Molecule オブジェクト → 相互作用検出（via interaction クラス群）
3. 検出結果 → IFP 格納
4. IFP → 解析出力（DataFrames, bitvectors, plots）

### Supported Input Formats

✅ **Supported** (Docking Analysis):
- PDB files: Protein structures（`Chem.MolFromPDBFile` 経由）
- SDF files: Ligand structures from docking（`sdf_supplier` クラス）
- MOL2 files: Limited support（`mol2_supplier` クラス）

⚠️ **Stub only** (raises NotImplementedError with migration guidance):
- PDBQT files: `pdbqt_supplier` クラスは存在するが `pdbqt_to_mol()` で NotImplementedError。PDB/SDF への変換が必要。
- MDAnalysis AtomGroup: `Molecule.from_mda()` は NotImplementedError

❌ **Not Supported** (MD Trajectory Analysis):
- Trajectory formats (XTC, TRR, DCD): MD シミュレーション解析
- `Fingerprint.run()` の実質的な軌跡解析: MDAnalysis 依存のため不可

### Key Dependencies

- **RDKit**: Chemical informatics and molecular representations（主要）
- **NetworkX**: Network analysis for interaction networks
- **Pandas**: Data manipulation and analysis
- **NumPy/SciPy**: Numerical computations
- **multiprocess/dill**: 並列処理

### Removed Dependencies

- **MDAnalysis**: ❌ Removed due to LGPL3 licensing. Functions requiring MDAnalysis now raise `NotImplementedError` with migration guidance.

## Pending Upstream Sync Tasks

upstream (chemosim-lab/ProLIF v2.1.0) との差分から、以下4タスクが未実装。  
詳細: `docs/upstream-sync-plan.md`

| # | タスク | 優先度 | 対象ファイル |
|---|--------|--------|------------|
| 1 | `numpy>=1.13.3,<2` → `numpy>=1.13.3`（NumPy 2.x 対応） | HIGH | `pyproject.toml` |
| 2 | `ResidueGroup.number` dtype: `uint16` → `uint32` | HIGH | `prolif/residue.py:217,221` |
| 3 | `split_molecule()` 関数の追加（RDKit のみ、MDAnalysis 不要） | HIGH | `prolif/molecule.py`, `prolif/__init__.py` |
| 4 | `get_n_jobs()` による並列数 cap 導入（`MAX_JOBS=10`） | MEDIUM | `prolif/parallel.py`, `prolif/fingerprint.py` |

完了後は `prolif/_version.py` を `2.2.0` に更新する。

## Migration Guide (MDAnalysis Removal)

### Hydrophobic Interaction Definition Changes

**⚠️ IMPORTANT CHANGE**: バージョン 2.1.0 で Hydrophobic 相互作用の定義を改良。  
N/O/F に結合した炭素原子を疎水性から除外するよう SMARTS を修正。

**Impact on Results:**
- Previous versions (≤2.0.3): Included aromatic carbons connected to N/O/F as hydrophobic
- Current version (≥2.1.0): Properly excludes all carbons connected to N/O/F from hydrophobic detection
- 結果として疎水性相互作用の数は減少するが科学的精度が向上

**Example:**
- Original ProLIF 2.0.3: 22 hydrophobic matches per ligand
- Current ProLIF 2.1.2: 14 hydrophobic matches per ligand (8 N/O/F-connected carbons excluded)

### Code Migration

**Old (MDAnalysis-based):**
```python
import MDAnalysis as mda
import prolif

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

protein_mol = Chem.MolFromPDBFile("protein.pdb", removeHs=False)
protein = prolif.Molecule.from_rdkit(protein_mol)

ligands = prolif.sdf_supplier("ligands.sdf")
fp = prolif.Fingerprint()
fp.run_from_iterable(ligands, protein)
```

### Stubs and Error Messages

MDAnalysis 依存の関数はすべて `NotImplementedError` を送出し、移行ガイダンスを提供する。

| Function | Status | Alternative |
|----------|--------|------------|
| `Molecule.from_mda()` | ⚠️ Stub (NotImplementedError) | Use `Molecule.from_rdkit()` |
| `pdbqt_supplier.pdbqt_to_mol()` | ⚠️ Stub (NotImplementedError) | Convert PDBQT to PDB/SDF first |
| `select_over_trajectory()` | ⚠️ Stub (NotImplementedError) | Use RDKit molecule methods |
| `Fingerprint.run()` (trajectory mode) | ⚠️ Partial stub | Use `run_from_iterable()` instead |
| Water bridge (trajectory mode) | ⚠️ Partial stub | Use `WaterBridge.run_from_iterable()` |

## Testing Structure

- Main tests in `tests/` directory
- Plotting tests in `tests/plotting/`
- Test fixtures defined in `conftest.py` with example molecules
- `tests/test_mdanalysis_removal.py`: MDAnalysis 除去後の NotImplementedError 動作を検証

## Data Files

- Example molecular data in `prolif/data/` (MOL2, PDB files)
- Vina docking outputs in `prolif/data/vina/`
- Network visualization assets in `prolif/plotting/network/`

## Code Style

- Uses ruff for linting and formatting
- Uses black for Jupyter notebook formatting
- Type hints enforced with mypy
- Follows PEP 8 style guidelines
