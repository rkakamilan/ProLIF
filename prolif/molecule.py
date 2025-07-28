"""
Reading proteins and ligands --- :mod:`prolif.molecule`
=======================================================
"""

import copy
from collections.abc import Iterable, Iterator, Sequence
from operator import attrgetter
from typing import TYPE_CHECKING, Any, TypeVar, Union, overload

try:
    import MDAnalysis as mda

    _HAS_MDANALYSIS = True
except ImportError:
    _HAS_MDANALYSIS = False

    # ダミーモジュールを定義
    class _DummyMDA:
        class SelectionError(Exception):
            pass

        def __getattr__(self, name):
            raise NotImplementedError("MDAnalysis is no longer supported")

    mda = _DummyMDA()

from rdkit import Chem

from prolif.rdkitmol import BaseRDKitMol
from prolif.residue import Residue, ResidueGroup
from prolif.utils import split_mol_by_residues

if TYPE_CHECKING:
    from pathlib import Path

    from prolif.typeshed import MDAObject, ResidueKey

    Self = TypeVar("Self", bound=BaseRDKitMol)


class Molecule(BaseRDKitMol):
    """Main molecule class that behaves like an RDKit :class:`~rdkit.Chem.rdchem.Mol`
    with extra attributes (see examples below). The main purpose of this class
    is to access residues as fragments of the molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        A ligand or protein with a single conformer
    use_segid: bool, default = False
        Use the segment number rather than the chain identifier as a chain.

    Attributes
    ----------
    residues : prolif.residue.ResidueGroup
        A dictionnary storing one/many :class:`~prolif.residue.Residue` indexed
        by :class:`~prolif.residue.ResidueId`. The residue list is sorted.
    n_residues : int
        Number of residues

    Examples
    --------

    .. ipython:: python
        :okwarning:

        import MDAnalysis as mda
        import prolif
        u = mda.Universe(prolif.datafiles.TOP, prolif.datafiles.TRAJ)
        mol = u.select_atoms("protein").convert_to("RDKIT")
        mol = prolif.Molecule(mol)
        mol

    You can also create a Molecule directly from a
    :class:`~MDAnalysis.core.universe.Universe`:

    .. ipython:: python
        :okwarning:

        mol = prolif.Molecule.from_mda(u, "protein")
        mol


    Notes
    -----
    Residues can be accessed easily in different ways:

    .. ipython:: python

        mol["TYR38.A"] # by resid string (residue name + number + chain)
        mol[42] # by index (from 0 to n_residues-1)
        mol[prolif.ResidueId("TYR", 38, "A")] # by ResidueId

    See :mod:`prolif.residue` for more information on residues

    .. versionchanged:: 2.1.0
        Added `use_segid`.
    """

    def __init__(self, mol: Chem.Mol, *, use_segid: bool = False) -> None:
        super().__init__(mol)
        # set mapping of atoms
        for atom in self.GetAtoms():
            atom.SetUnsignedProp("mapindex", atom.GetIdx())
        # split in residues
        residues = split_mol_by_residues(self)
        residues = [Residue(mol, use_segid=use_segid) for mol in residues]
        residues.sort(key=attrgetter("resid"))
        self.residues = ResidueGroup(residues)

    @classmethod
    def from_mda(
        cls,
        obj: "MDAObject",
        selection: str | None = None,
        *,
        use_segid: bool | None = None,
        **kwargs: Any,
    ) -> "Molecule":
        """Creates a Molecule from an MDAnalysis object (DEPRECATED)

        .. deprecated::
            MDAnalysis support has been removed. Use alternative methods:
            - For PDB files: use :meth:`from_file` or RDKit directly
            - For single structures: use supplier classes (sdf_supplier, pdbqt_supplier)
            - See migration guide in documentation for details.

        Raises
        ------
        NotImplementedError
            This method is no longer supported
        """
        raise NotImplementedError(
            "Molecule.from_mda() is no longer supported due to removal of "
            "MDAnalysis dependency. Alternative approaches:\n"
            "  - For PDB files: use "
            "Molecule.from_rdkit(Chem.MolFromPDBFile('file.pdb'))\n"
            "  - For SDF files: use sdf_supplier('file.sdf')\n"
            "  - For PDBQT files: use pdbqt_supplier(['file.pdbqt'], template)\n"
            "  - For single structures: use RDKit directly\n"
            "See documentation for detailed migration guide."
        )

    @classmethod
    def _use_segid(cls, obj: "MDAObject", use_segid: bool | None = None) -> bool:
        """Whether to use the segment index or the chainID as a chain (DEPRECATED)."""
        raise NotImplementedError(
            "_use_segid() is no longer supported due to removal of "
            "MDAnalysis dependency."
        )

    @classmethod
    def from_rdkit(
        cls,
        mol: Chem.Mol,
        resname: str = "UNL",
        resnumber: int = 1,
        chain: str = "",
        use_segid: bool = False,
    ) -> "Molecule":
        """Creates a Molecule from an RDKit molecule

        While directly instantiating a molecule with ``prolif.Molecule(mol)``
        would also work, this method insures that every atom is linked to an
        AtomPDBResidueInfo which is required by ProLIF

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            The input RDKit molecule
        resname : str
            The default residue name that is used if none was found
        resnumber : int
            The default residue number that is used if none was found
        chain : str
            The default chain Id that is used if none was found
        use_segid: bool, default = False
            Use the segment number rather than the chain identifier as a chain

        Notes
        -----
        This method only checks for an existing AtomPDBResidueInfo in the first
        atom. If none was found, it will patch all atoms with the one created
        from the method's arguments (resname, resnumber, chain).

        .. versionchanged:: 2.1.0
            Added `use_segid`.
        """
        if mol.GetAtomWithIdx(0).GetMonomerInfo():
            return cls(mol, use_segid=use_segid)
        mol = copy.deepcopy(mol)
        for atom in mol.GetAtoms():
            mi = Chem.AtomPDBResidueInfo(
                f" {atom.GetSymbol():<3.3}",
                residueName=resname,
                residueNumber=resnumber,
                chainId=chain,
            )
            atom.SetMonomerInfo(mi)
        return cls(mol)

    def __iter__(self) -> Iterator[Residue]:
        yield from self.residues.values()

    def __getitem__(self, key: "ResidueKey") -> Residue:
        return self.residues[key]

    def __repr__(self) -> str:  # pragma: no cover
        name = ".".join([self.__class__.__module__, self.__class__.__name__])
        params = f"{self.n_residues} residues and {self.GetNumAtoms()} atoms"
        return f"<{name} with {params} at {id(self):#x}>"

    @property
    def n_residues(self) -> int:
        return len(self.residues)


class pdbqt_supplier(Sequence[Molecule]):
    """Supplies molecules, given paths to PDBQT files

    Parameters
    ----------
    paths : list
        A list (or any iterable) of PDBQT files
    template : rdkit.Chem.rdchem.Mol
        A template molecule with the correct bond orders and charges. It must
        match exactly the molecule inside the PDBQT file.
    converter_kwargs : dict | None
        Keyword arguments passed to the RDKitConverter of MDAnalysis
    resname : str
        Residue name for every ligand
    resnumber : int
        Residue number for every ligand
    chain : str
        Chain ID for every ligand


    Returns
    -------
    suppl : Sequence
        A sequence that provides :class:`Molecule` objects

    Example
    -------
    The supplier is typically used like this::

        >>> import glob
        >>> pdbqts = glob.glob("docking/ligand1/*.pdbqt")
        >>> lig_suppl = pdbqt_supplier(pdbqts, template)
        >>> for lig in lig_suppl:
        ...     # do something with each ligand

    .. versionchanged:: 1.0.0
        Molecule suppliers are now sequences that can be reused, indexed,
        and can return their length, instead of single-use generators.

    .. versionchanged:: 1.1.0
        Because the PDBQT supplier needs to strip hydrogen atoms before
        assigning bond orders from the template, it used to replace them
        entirely with hydrogens containing new coordinates. It now directly
        uses the hydrogen atoms present in the file and won't add explicit
        ones anymore, to prevent the fingerprint from detecting hydrogen bonds
        with "random" hydrogen atoms.
        A lot of irrelevant warnings and logs have been disabled as well.

    """

    def __init__(
        self,
        paths: Iterable[Union[str, "Path"]],
        template: Chem.Mol,
        converter_kwargs: dict | None = None,
        **kwargs: Any,
    ) -> None:
        self.paths = list(paths)
        self.template = template
        converter_kwargs = converter_kwargs or {}
        converter_kwargs.pop("NoImplicit", None)
        self.converter_kwargs = converter_kwargs
        self._kwargs = kwargs

    def __iter__(self) -> Iterator[Molecule]:
        # Check if MDAnalysis is available at iteration time
        try:
            import MDAnalysis
        except ImportError:
            raise NotImplementedError(
                "PDBQT file support requires MDAnalysis dependency, which has been "
                "removed. Alternative approaches:\n"
                "  - Convert PDBQT files to PDB/SDF format using external tools\n"
                "  - Use sdf_supplier() for SDF ligand files\n"
                "  - Use RDKit directly: Molecule.from_rdkit(Chem.MolFromPDBFile())\n"
                "See documentation for detailed migration guide."
            ) from None

        for pdbqt_path in self.paths:
            yield self.pdbqt_to_mol(pdbqt_path)

    @overload
    def __getitem__(self, index: int) -> Molecule: ...
    @overload
    def __getitem__(self, index: slice) -> "pdbqt_supplier": ...
    def __getitem__(self, index: int | slice) -> Union[Molecule, "pdbqt_supplier"]:
        if isinstance(index, slice):
            return pdbqt_supplier(
                self.paths[index],
                self.template,
                converter_kwargs=self.converter_kwargs,
                **self._kwargs,
            )

        pdbqt_path = self.paths[index]
        return self.pdbqt_to_mol(pdbqt_path)

    def pdbqt_to_mol(self, pdbqt_path: Union[str, "Path"]) -> Molecule:
        raise NotImplementedError(
            "PDBQT file support requires MDAnalysis, which is no longer available. "
            "Please convert PDBQT files to PDB or SDF format using external tools "
            "(e.g., Open Babel: obabel -ipdbqt file.pdbqt -opdb -O file.pdb)."
        )

    @staticmethod
    def _adjust_hydrogens(template: Chem.Mol, pdbqt_mol: Chem.Mol) -> Chem.Mol:
        """Adjust hydrogens (DEPRECATED - MDAnalysis functionality removed)"""
        raise NotImplementedError(
            "_adjust_hydrogens is no longer supported due to removal of "
            "MDAnalysis dependency. PDBQT functionality has been removed."
        )

    def __len__(self) -> int:
        return len(self.paths)


class sdf_supplier(Sequence[Molecule]):
    """Supplies molecules, given a path to an SDFile

    Parameters
    ----------
    path : str
        A path to the .sdf file
    sanitize : bool
        Whether to sanitize each molecule or not.
    resname : str
        Residue name for every ligand
    resnumber : int
        Residue number for every ligand
    chain : str
        Chain ID for every ligand

    Returns
    -------
    suppl : Sequence
        A sequence that provides :class:`Molecule` objects. Can be indexed

    Example
    -------
    The supplier is typically used like this::

        >>> lig_suppl = sdf_supplier("docking/output.sdf")
        >>> for lig in lig_suppl:
        ...     # do something with each ligand

    .. versionchanged:: 1.0.0
        Molecule suppliers are now sequences that can be reused, indexed,
        and can return their length, instead of single-use generators.

    .. versionchanged:: 2.1.0
        Added ``sanitize`` parameter (defaults to ``True``, same behavior as before).

    """

    def __init__(self, path: str, sanitize: bool = True, **kwargs: Any) -> None:
        self.path = path
        self._suppl = Chem.SDMolSupplier(path, removeHs=False, sanitize=sanitize)
        self._kwargs = kwargs

    def __iter__(self) -> Iterator[Molecule]:
        for mol in self._suppl:
            yield Molecule.from_rdkit(mol, **self._kwargs)

    @overload
    def __getitem__(self, index: int) -> Molecule: ...
    @overload
    def __getitem__(self, index: slice) -> "sdf_supplier": ...
    def __getitem__(self, index: int | slice) -> Union[Molecule, "sdf_supplier"]:
        if isinstance(index, slice):
            suppl = sdf_supplier(self.path, **self._kwargs)
            suppl._suppl = [  # type: ignore[assignment]
                self[i] for i in range(*index.indices(len(self)))
            ]
            return suppl

        mol = self._suppl[index]
        return Molecule.from_rdkit(mol, **self._kwargs)

    def __len__(self) -> int:
        return len(self._suppl)


class mol2_supplier(Sequence[Molecule]):
    """Supplies molecules, given a path to a MOL2 file

    Parameters
    ----------
    path : str
        A path to the .mol2 file
    sanitize : bool
        Whether to sanitize each molecule or not.
    cleanup_substructures : bool
        Toggles standardizing some substructures found in mol2 files, based on atom
        types.
    resname : str
        Residue name for every ligand
    resnumber : int
        Residue number for every ligand
    chain : str
        Chain ID for every ligand

    Returns
    -------
    suppl : Sequence
        A sequence that provides :class:`Molecule` objects

    Example
    -------
    The supplier is typically used like this::

        >>> lig_suppl = mol2_supplier("docking/output.mol2")
        >>> for lig in lig_suppl:
        ...     # do something with each ligand

    .. versionchanged:: 1.0.0
        Molecule suppliers are now sequences that can be reused, indexed,
        and can return their length, instead of single-use generators.

    .. versionchanged:: 2.1.0
        Added ``cleanup_substructures`` and ``sanitize`` parameters
        (default to ``True``, same behavior as before).

    """

    def __init__(
        self,
        path: Union[str, "Path"],
        cleanup_substructures: bool = True,
        sanitize: bool = True,
        **kwargs: Any,
    ) -> None:
        self.path = path
        self.cleanup_substructures = cleanup_substructures
        self.sanitize = sanitize
        self._kwargs = kwargs

    def __iter__(self) -> Iterator[Molecule]:
        block: list[str] = []
        with open(self.path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                if block and line.startswith("@<TRIPOS>MOLECULE"):
                    yield self.block_to_mol(block)
                    block = []
                block.append(line)
            yield self.block_to_mol(block)

    def block_to_mol(self, block: list[str]) -> Molecule:
        mol = Chem.MolFromMol2Block(
            "".join(block),
            removeHs=False,
            cleanupSubstructures=self.cleanup_substructures,
            sanitize=self.sanitize,
        )
        return Molecule.from_rdkit(mol, **self._kwargs)

    @overload
    def __getitem__(self, index: int) -> Molecule: ...
    @overload
    def __getitem__(self, index: slice) -> "mol2_supplier": ...
    def __getitem__(self, index: int | slice) -> Union[Molecule, "mol2_supplier"]:
        if isinstance(index, slice):
            raise NotImplementedError("Slicing not available for mol2_supplier.")

        if index < 0:
            index %= len(self)
        mol_index = -1
        molblock_started = False
        block: list[str] = []
        with open(self.path) as f:
            for line in f:
                if line.startswith("@<TRIPOS>MOLECULE"):
                    mol_index += 1
                    if mol_index > index:
                        return self.block_to_mol(block)
                    if mol_index == index:
                        molblock_started = True
                if molblock_started:
                    block.append(line)
        if block:
            return self.block_to_mol(block)
        raise ValueError(f"Could not parse molecule with index {index}")

    def __len__(self) -> int:
        with open(self.path) as f:
            return sum(line.startswith("@<TRIPOS>MOLECULE") for line in f)
