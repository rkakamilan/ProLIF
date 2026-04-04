# ProLIF MDAnalysis除去設計書 (最小改修版)

## 背景・目的

ProLIFからMDAnalysisの依存性を除去し、ドッキング結果解析に特化する。新しいクラス実装は避け、既存コードの最小限の修正で対応する。

## 方針

### 基本戦略
1. **MDAnalysis依存部分のエラー対応**: トラジェクトリ関連機能で適切なエラーメッセージを表示
2. **VdW半径の独自実装**: MDAnalysisのvdwradiiデータを独自実装に置き換え
3. **既存APIの維持**: 可能な限り既存のAPIを保持
4. **段階的移行**: deprecation warningで段階的に移行

### 対象機能の分類

#### 完全に廃止（エラー対応）
- `Molecule.from_mda()` → ImportErrorまたはNotImplementedError
- `Fingerprint.run()` with trajectory → エラーメッセージ
- `select_over_trajectory()` → エラーメッセージ
- `UpdatingAtomGroup`関連 → 代替実装またはエラー

#### 独自実装で置き換え
- VdW半径データ (`prolif/interactions/constants.py`)
- 基本的な原子データ

#### そのまま維持
- RDKitベースの単一構造処理
- 既存の相互作用検出ロジック
- `sdf_supplier`, `pdbqt_supplier`など

## 具体的な実装方針

### 1. VdW半径データの独自実装

#### 現在のコード (`prolif/interactions/constants.py`)
```python
from MDAnalysis.topology.tables import vdwradii
```

#### 修正後
```python
# 文献ベースのVdW半径データを独自定義
VDW_RADII = {
    'H': 1.2,   # van der Waals radius in Angstroms
    'He': 1.4,
    'Li': 1.82,
    'Be': 1.53,
    'B': 1.92,
    'C': 1.7,
    'N': 1.55,
    'O': 1.52,
    'F': 1.47,
    'Ne': 1.54,
    # ... 完全なデータセット
}

# MDAnalysisのvdwradiiと同じインターフェースを提供
vdwradii = VDW_RADII  # 後方互換性のため
```

### 2. MDAnalysis依存関数のエラー対応

#### `prolif/molecule.py`の修正

```python
@classmethod
def from_mda(cls, obj, selection=None, **kwargs):
    """MDAnalysisオブジェクトからMoleculeを作成（廃止予定）"""
    raise NotImplementedError(
        "Molecule.from_mda() is no longer supported. "
        "Please use from_file() for PDB files, or the supplier classes "
        "for SDF/MOL2/PDBQT files. "
        "See documentation for migration guide."
    )
```

#### `prolif/fingerprint.py`の修正

```python
def run(self, trajectory, ligand_mol, protein_mol, **kwargs):
    """フィンガープリント計算を実行"""
    # トラジェクトリかどうかを判定
    if hasattr(trajectory, '__iter__') and not isinstance(trajectory, (str, Path)):
        # トラジェクトリの場合はエラー
        raise NotImplementedError(
            "Trajectory analysis is no longer supported. "
            "Please use analyze_structures() for single structure analysis. "
            "For multiple structures, use the supplier classes."
        )
    
    # 既存の単一構造解析コードはそのまま
    return self._analyze_single_structure(ligand_mol, protein_mol, **kwargs)
```

#### `prolif/utils.py`の修正

```python
def select_over_trajectory(*args, **kwargs):
    """トラジェクトリ上での原子選択（廃止）"""
    raise NotImplementedError(
        "select_over_trajectory() is no longer supported. "
        "This function was designed for MD trajectory analysis."
    )
```

### 3. 水分子ブリッジの簡略化

#### `prolif/interactions/water_bridge.py`の修正

```python
try:
    from MDAnalysis.core.groups import UpdatingAtomGroup
    _HAS_MDANALYSIS = True
except ImportError:
    _HAS_MDANALYSIS = False
    
    # ダミークラスを定義
    class UpdatingAtomGroup:
        def __init__(self, *args, **kwargs):
            raise NotImplementedError(
                "UpdatingAtomGroup requires MDAnalysis, which is no longer supported. "
                "Water bridge analysis for trajectories is not available."
            )

class WaterBridge:
    def __init__(self, **kwargs):
        if not _HAS_MDANALYSIS and 'updating' in kwargs:
            raise NotImplementedError(
                "Dynamic water bridge analysis requires MDAnalysis, "
                "which is no longer supported."
            )
        # 静的解析のみサポート
```

### 4. パッケージレベルでの依存関係除去

#### `pyproject.toml`の修正

```toml
dependencies = [
    "pandas>=1.1.0",
    "numpy>=1.13.3,<2",
    "scipy>=1.3.0",
    # "mdanalysis>=2.2.0,<3; python_version<'3.13'",  # 削除
    # "mdanalysis>=2.7.0,<3; python_version>='3.13'",  # 削除
    "networkx>=2.5.0",
    "tqdm",
    "multiprocess",
    "dill",
]
```

### 5. 推奨される新しい使用パターン

#### 既存コード（廃止予定）
```python
import MDAnalysis as mda
import prolif

u = mda.Universe("protein.pdb", "trajectory.xtc")
ligand = u.select_atoms("resname LIG")
protein = u.select_atoms("protein")

fp = prolif.Fingerprint()
fp.run(u.trajectory, ligand, protein)
```

#### 新しいコード（推奨）
```python
import prolif

# 単一構造解析
protein = prolif.Molecule.from_file("protein.pdb")
ligand = prolif.Molecule.from_file("ligand.sdf")

fp = prolif.Fingerprint()
ifp = fp.analyze_structures(ligand, protein)

# 複数構造（ドッキング結果）
ligand_suppl = prolif.sdf_supplier("docking_results.sdf")
results = []
for ligand in ligand_suppl:
    ifp = fp.analyze_structures(ligand, protein)
    results.append(ifp)
```

## 移行計画

### フェーズ1: 基本的な置き換え（1週間）
1. VdW半径データの独自実装
2. MDAnalysis import部分の条件分岐追加
3. 基本的なエラーメッセージの実装

### フェーズ2: 関数レベルの対応（1週間）
1. `from_mda()`のエラー実装
2. `run()`メソッドのトラジェクトリ判定追加
3. `select_over_trajectory()`のエラー実装

### フェーズ3: テストと検証（1週間）
1. 既存テストの修正
2. エラーメッセージの確認
3. ドッキング解析での動作確認

## 利点

1. **最小限の変更**: 既存コードの大部分をそのまま維持
2. **明確なエラー**: 利用不可機能で適切なエラーメッセージ
3. **段階的移行**: 既存ユーザーへの影響を最小化
4. **保守性**: 複雑な新クラス実装を避けることで保守負荷軽減

## ファイル別の具体的な修正箇所

### 必須修正ファイル
1. `pyproject.toml` - MDAnalysis依存関係の除去
2. `prolif/interactions/constants.py` - VdW半径の独自実装
3. `prolif/molecule.py` - `from_mda()`のエラー実装
4. `prolif/fingerprint.py` - `run()`のトラジェクトリ判定
5. `prolif/utils.py` - `select_over_trajectory()`のエラー実装
6. `prolif/interactions/water_bridge.py` - 条件分岐の追加

### テスト修正
1. MDAnalysisを使用するテストの無効化またはスキップ
2. エラーメッセージのテスト追加
3. 単一構造解析テストの強化

この最小改修アプローチにより、コードベースへの影響を最小限に抑えながら、MDAnalysis依存関係を効果的に除去できます。