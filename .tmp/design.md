# ProLIF MDAnalysis除去設計書

## 背景・目的

ProLIFは現在LGPL3ライセンスのMDAnalysisに依存しており、これがより制限の少ないライセンス（Apache 2.0など）での利用を妨げている。本設計では、MDAnalysisの依存性を完全に除去し、代替ソリューションに移行する方法を定義する。

## 現在のMDAnalysis使用状況分析

### 主要使用箇所

1. **分子ファイル読み込み** (`prolif/molecule.py`)
   - PDB、XTC、PDBQTファイルの読み込み
   - `ag.convert_to.rdkit()`でRDKitMolへの変換
   - トポロジー情報の取得

2. **トラジェクトリ処理** (`prolif/fingerprint.py`, `prolif/utils.py`)
   - `MDAnalysis.AtomGroup`、`Universe`、`Timestep`の使用
   - トラジェクトリイテレーション
   - 原子選択とフィルタリング

3. **座標変換** (`prolif/interactions/constants.py`)
   - van der Waals半径の取得
   - `MDAnalysis.topology.tables.vdwradii`

4. **水分子ブリッジ** (`prolif/interactions/water_bridge.py`)
   - `UpdatingAtomGroup`の使用

5. **ユーティリティ関数** (`prolif/utils.py`)
   - 原子選択とトラジェクトリ処理
   - `select_over_trajectory`関数

## 代替ライブラリ選定

### MDTraj（推奨）

**利点:**
- BSD-3-Clauseライセンス（制限が少ない）
- 豊富なファイル形式対応（PDB、XTC、DCD、NetCDF、HDF5など）
- NumPy/SciPyベースで高速
- RDKitとの統合例も存在
- 活発な開発とメンテナンス

**対応可能機能:**
- トラジェクトリ読み込み・処理
- 座標データアクセス
- トポロジー情報の取得
- 原子選択（制限あり）

### 追加検討ライブラリ

1. **独自実装 + RDKit**
   - PDB読み込みはRDKitで対応可能
   - 単純なファイル形式は独自パーサーで実装

2. **PDBファイル用途特化**
   - BioPythonのPDB.PDBParser（Biopython License）
   - 軽量な独自PDBパーサー

## アーキテクチャ設計案

### 新しいモジュール構造

```
prolif/
├── io/                    # 新設：ファイルI/Oモジュール
│   ├── __init__.py
│   ├── base.py           # 基底クラス
│   ├── mdtraj_io.py      # MDTrajベースのI/O
│   ├── rdkit_io.py       # RDKitベースのI/O
│   └── parsers.py        # 独自パーサー（必要に応じて）
├── trajectory/           # 新設：トラジェクトリ処理
│   ├── __init__.py
│   ├── base.py          # 基底クラス
│   ├── mdtraj_traj.py   # MDTrajベースの実装
│   └── utils.py         # トラジェクトリユーティリティ
└── constants/           # 新設：定数とデータ
    ├── __init__.py
    ├── vdw_radii.py     # VdW半径データ
    └── atom_data.py     # 原子データ
```

### 抽象化レイヤー設計

#### 1. TrajectoryReader基底クラス

```python
from abc import ABC, abstractmethod
from typing import Iterator, Optional, Union
import numpy as np

class TrajectoryReader(ABC):
    @abstractmethod
    def __init__(self, filename: str, topology: Optional[str] = None):
        pass
    
    @abstractmethod
    def __iter__(self) -> Iterator['Frame']:
        pass
    
    @abstractmethod
    def __len__(self) -> int:
        pass
    
    @abstractmethod
    def __getitem__(self, index: Union[int, slice]) -> 'Frame':
        pass
    
    @property
    @abstractmethod
    def n_atoms(self) -> int:
        pass
    
    @property
    @abstractmethod
    def n_frames(self) -> int:
        pass
```

#### 2. Frame データクラス

```python
@dataclass
class Frame:
    coordinates: np.ndarray  # shape: (n_atoms, 3)
    time: Optional[float] = None
    box: Optional[np.ndarray] = None  # unit cell
    topology: Optional['Topology'] = None
```

#### 3. Topology クラス

```python
class Topology:
    def __init__(self, rdkit_mol: Chem.Mol):
        self._mol = rdkit_mol
        # MDAnalysisのトポロジー情報をRDKitベースで再構築
    
    def select_atoms(self, selection: str) -> np.ndarray:
        # 基本的な原子選択機能
        pass
    
    @property
    def atoms(self) -> List['Atom']:
        pass
```

### 移行戦略

#### フェーズ1: I/Oレイヤーの実装
1. MDTrajベースのTrajectoryReaderクラス実装
2. RDKitベースの単一構造読み込み実装
3. 抽象化インターフェースの定義

#### フェーズ2: 既存コードの移行
1. `prolif.molecule.Molecule.from_mda()` → `from_file()`
2. `prolif.fingerprint.Fingerprint.run()`のトラジェクトリ処理部分
3. `prolif.utils`のユーティリティ関数

#### フェーズ3: 高度な機能の実装
1. 水分子ブリッジ処理の再実装
2. 原子選択機能の拡張
3. パフォーマンス最適化

## 互換性戦略

### API互換性の維持

既存のユーザーAPIを可能な限り維持:

```python
# 現在のAPI
mol = prolif.Molecule.from_mda(u, "protein")

# 新しいAPI（後方互換性あり）
mol = prolif.Molecule.from_file("protein.pdb")
mol = prolif.Molecule.from_trajectory(traj_reader, frame_index=0)

# 内部的にはMDTrajを使用
```

### 段階的移行

1. **deprecation warning**の追加
2. 両方のバックエンドをサポート（optional dependency）
3. 最終的にMDAnalysisサポートを削除

## リスク分析

### 高リスク項目

1. **原子選択機能の制限**
   - MDAnalysisの高度な選択構文が使用不可
   - 対策: 基本的な選択機能を独自実装

2. **トラジェクトリ形式の対応不足**
   - MDTraj未対応の形式が存在する可能性
   - 対策: 段階的に対応形式を拡張

3. **パフォーマンスの劣化**
   - I/Oボトルネックの可能性
   - 対策: ベンチマークとプロファイリング

### 中リスク項目

1. **メモリ使用量の変化**
   - MDTrajとMDAnalysisでメモリ効率が異なる
   - 対策: メモリ使用量の継続監視

2. **依存関係の増加**
   - MDTraj、その他のライブラリ追加
   - 対策: 最小限の依存関係に留める

## 実装優先度

### 高優先度
1. 基本的なPDB/XTC読み込み機能
2. `Molecule.from_mda()`の代替実装
3. 単一フレーム処理

### 中優先度
1. マルチフレームトラジェクトリ処理
2. 原子選択機能
3. パフォーマンス最適化

### 低優先度
1. 高度なトラジェクトリ分析機能
2. MDAnalysis特有の機能の代替実装
3. 追加ファイル形式のサポート

## 成功指標

1. **機能カバレッジ**: 既存機能の95%以上を維持
2. **パフォーマンス**: 既存実装と同等以上の性能
3. **API互換性**: 既存ユーザーコードの90%以上が無修正で動作
4. **テストカバレッジ**: 全テストの95%以上がパス
5. **ライセンス**: LGPL3依存関係の完全除去