# Upstream Sync Plan: chemosim-lab/ProLIF → prolif-lite

**作成日**: 2026-04-04  
**対象upstream**: https://github.com/chemosim-lab/ProLIF (v2.1.0, 2026-01-20)  
**当リポジトリ**: prolif-lite v2.1.3 (MDAnalysis除去フォーク)

---

## 調査サマリー

upstream v2.1.0 で追加された機能のうち、**MDAnalysis不要のもの**はすべてバックポート対象。MDAnalysis依存の機能（TrajectoryPoolQueue, parallel_strategy など）は除外。

---

## 差分一覧

### ✅ 既に実装済み（対応不要）

| 機能 | upstream PR/Issue |
|------|------------------|
| Hydrophobic SMARTS修正（N/O/F接続炭素の除外） | v2.1.0 |
| HBond SMARTS修正（帯電芳香族窒素等の除外） | v2.1.0 |
| VdWContact `preset` パラメータ（mdanalysis/rdkit/csd） | v2.1.0 |
| IFP.interactions() イテレータ | v2.1.0 |
| BridgedInteraction / WaterBridge（RDKit版） | v2.1.0 |
| Interaction.all(), any(), best() メソッド | v2.1.0 |
| Complex3D: save_png, only_interacting, remove_hydrogens | v2.1.0 |
| LigNetwork: show_interaction_data, save_png | v2.1.0 |
| plot_barcode() メソッド | v2.1.0 |
| ResidueId: 空白除去、セグメントID対応 | v2.1.0 |
| sdf_supplier: sanitize パラメータ | v2.1.0 |
| use_segid パラメータ（Molecule, Fingerprint） | v2.1.0 |

### ❌ 未実装・要対応

以下の4項目が未実装または不正確な状態にある。

---

## 実装タスク（優先度順）

---

### Task 1: numpy上限制約の撤廃【優先度: HIGH / 工数: XS】

**問題**:  
`pyproject.toml` に `numpy>=1.13.3,<2` とあり、NumPy 2.x を人工的にブロックしている。  
upstream は `numpy>=1.13.3`（上限なし）に変更済み。NumPy 2.x はすでに広く使われており、この制約はユーザーの環境構築を不必要に制限する。

**変更ファイル**: `pyproject.toml`

**変更内容**:
```diff
- "numpy>=1.13.3,<2",
+ "numpy>=1.13.3",
```

**リスク**: 低。NumPy 2.x の破壊的変更はほぼないが、念のためテスト通過を確認する。

---

### Task 2: ResidueGroup の number 型を uint16 → uint32 に修正【優先度: HIGH / 工数: XS】

**問題**:  
`prolif/residue.py` 内の `ResidueGroup` クラスで `number` 配列の dtype が `uint16`（最大65535）になっている。  
upstream はこれを `uint32` に修正済み。大きなシステム（分子動力学シミュレーション等）では残基番号が65535を超える場合があり、`uint16` はオーバーフローを引き起こす。

**変更ファイル**: `prolif/residue.py`

**現状** (line 217-221):
```python
self.number = np.array([], dtype=np.uint8)  # 空の場合
...
self.number = np.asarray(number, dtype=np.uint16)  # 非空の場合
```

**変更内容**:
```diff
- self.number = np.array([], dtype=np.uint8)
+ self.number = np.array([], dtype=np.uint32)
  ...
- self.number = np.asarray(number, dtype=np.uint16)
+ self.number = np.asarray(number, dtype=np.uint32)
```

**リスク**: 低。メモリ使用量がわずかに増えるが機能的影響なし。

---

### Task 3: `split_molecule` 関数の追加【優先度: HIGH / 工数: S】

**問題**:  
upstream v2.1.0 で追加された `split_molecule()` 関数が当リポジトリに存在しない。  
この関数はMDAnalysis不要で、RDKitのみで動作する。溶媒和複合体からリガンドや水分子を分離する実用的なユーティリティ。

**変更ファイル**: 
- `prolif/molecule.py`（関数実装）
- `prolif/__init__.py`（エクスポート追加）

**実装内容** (`prolif/molecule.py` 末尾に追加):
```python
def split_molecule(
    mol: "Molecule", predicate: "Callable[[ResidueId], bool]"
) -> "tuple[Molecule, Molecule]":
    """Splits a molecule into two based on a predicate function.

    The first molecule returned contains all residues of the input mol for which
    the predicate function was true, the second molecule contains the rest.
    This function is typically used to extract a ligand or water molecules from a file
    containing a solvated complex::

        >>> solvated_system = Molecule(...)
        >>> water_mol, protein_mol = split_molecule(
        ...     solvated_system, lambda x: x.name == "WAT"
        ... )

    Parameters
    ----------
    mol : prolif.molecule.Molecule
        Input molecule containing multiple residues
    predicate : Callable[[ResidueId], bool]
        A function that takes a ResidueId and returns True for residues to
        include in the first output molecule.

    Returns
    -------
    tuple[Molecule, Molecule]
        Two molecules: the first contains residues where predicate returned True,
        the second contains the remaining residues.

    .. versionadded:: 2.1.0

    """
    from rdkit import Chem

    lhs_atoms = []
    rhs_atoms = []
    for residue in mol:
        if predicate(residue.resid):
            lhs_atoms.extend(
                atom.GetUnsignedProp("mapindex") for atom in residue.GetAtoms()
            )
        else:
            rhs_atoms.extend(
                atom.GetUnsignedProp("mapindex") for atom in residue.GetAtoms()
            )

    def _extract(mol: "Molecule", atom_indices: list[int]) -> "Molecule":
        rw = Chem.RWMol(mol)
        # Remove atoms not in the list, in reverse order to keep indices valid
        remove = sorted(
            set(range(mol.GetNumAtoms())) - set(atom_indices), reverse=True
        )
        for idx in remove:
            rw.RemoveAtom(idx)
        return Molecule(rw.GetMol())

    return _extract(mol, lhs_atoms), _extract(mol, rhs_atoms)
```

**`prolif/__init__.py` の変更**:
```diff
 from prolif.molecule import Molecule, mol2_supplier, pdbqt_supplier, sdf_supplier
+from prolif.molecule import split_molecule
 from prolif.utils import (
```

**注意**: upstream の実装は `RWMol` の `RemoveAtom` を `mapindex` ベースで行っているが、インデックスが変動するため実装に注意が必要。upstream コードを精査して正確に移植する。

**リスク**: 中。`mapindex` の扱いを誤るとサイレントなバグになる可能性がある。テスト必須。

---

### Task 4: `get_n_jobs` 関数による並列数 cap の導入【優先度: MEDIUM / 工数: S】

**問題**:  
upstream v2.1.0 で `get_n_jobs(n_jobs, capped=True/False)` 関数が追加され、並列処理数の上限を `MAX_JOBS=10` で制限するようになった。  

**現状（当フォーク）**:  
`prolif/fingerprint.py` の `run()`（line 652-653）と `run_from_iterable()`（line 871-872）は `PROLIF_N_JOBS` 環境変数の読み込みを既に実装済み。**不足しているのは `MAX_JOBS=10` の cap 機構のみ**。大規模なドッキング解析（数千ポーズ等）で `n_jobs=None` を指定すると全CPUを使い切り、メモリ枯渇が起きる可能性がある。

**変更ファイル**:
- `prolif/parallel.py`（`get_n_jobs` 関数を新規追加）
- `prolif/fingerprint.py`（`n_jobs` 解決ロジックを `get_n_jobs` に委譲）

> `psutil` は upstream では追加されているが、Python 標準の `multiprocessing.cpu_count()` で代替可能。依存追加は不要。

**`prolif/parallel.py` に追加する実装**:
```python
import os
from multiprocessing import cpu_count

MAX_JOBS = 10

def get_n_jobs(n_jobs: int | None, *, capped: bool = True) -> int | None:
    """Determine the number of parallel jobs to use.

    Priority: explicit n_jobs > PROLIF_N_JOBS env var > logical CPU count.
    When capped=True, limits to MAX_JOBS to prevent excessive memory usage.

    Parameters
    ----------
    n_jobs : int or None
        Explicitly requested number of jobs. None means use all available.
    capped : bool
        Whether to cap at MAX_JOBS (default True for run(), False for run_from_iterable()).

    Returns
    -------
    int or None
        Number of jobs to use, or None if sequential.
    """
    if n_jobs is not None:
        return n_jobs
    env_jobs = os.environ.get("PROLIF_N_JOBS")
    if env_jobs is not None:
        return int(env_jobs) or None
    n_cpus = cpu_count()
    if capped:
        return min(n_cpus, MAX_JOBS)
    return n_cpus
```

**`prolif/fingerprint.py` の変更**:  
`run()` の lines 652-653 と `run_from_iterable()` の lines 871-872 の inline `n_jobs` 解決ロジックを `get_n_jobs()` 呼び出しに置き換える。

```diff
# run() 内 (line 652)
- if n_jobs is None:
-     n_jobs = int(os.environ.get("PROLIF_N_JOBS", "0")) or None
+ n_jobs = get_n_jobs(n_jobs, capped=True)

# run_from_iterable() 内 (line 871)
- if n_jobs is None:
-     n_jobs = int(os.environ.get("PROLIF_N_JOBS", "0")) or None
+ n_jobs = get_n_jobs(n_jobs, capped=False)
```

**リスク**: 低。後方互換性あり（明示的に `n_jobs` を指定している場合はそのまま使われる）。

---

## 実装の順序と推奨アプローチ

```
Phase 1 (即時・低リスク):
  ├── Task 1: numpy上限撤廃     (pyproject.toml 1行変更)
  └── Task 2: uint16 → uint32   (residue.py 2行変更)

Phase 2 (機能追加):
  ├── Task 3: split_molecule     (molecule.py + __init__.py)
  └── Task 4: get_n_jobs         (parallel.py + fingerprint.py + pyproject.toml)
```

各タスクは独立したブランチ・PRで対応する。

---

## スコープ外（対応しない）

| 機能 | 理由 |
|------|------|
| TrajectoryPoolQueue | MDAnalysis依存。当リポジトリの方針と矛盾 |
| parallel_strategy パラメータ | TrajectoryPoolQueueに依存 |
| select_over_trajectory の完全実装 | MDAnalysis依存。現状のスタブ（NotImplementedError）を維持 |
| WaterBridge.run() のtrajectory対応 | MDAnalysis依存 |

---

## テスト戦略

各タスクに対して以下のテストを追加・確認する:

- **Task 1**: `poe test` で既存テスト全通過を確認
- **Task 2**: 残基番号が65535を超えるケースのユニットテスト追加
- **Task 3**: `split_molecule` の単体テスト（正常系・空の結果・全選択）
- **Task 4**: `get_n_jobs` の単体テスト（環境変数, None, 明示値, capping）

---

## バージョン管理

実装完了後、`prolif/_version.py` を `2.2.0` に更新する。  
（機能追加を含むため、マイナーバージョンを上げる）

```python
__version__ = "2.2.0"
```
