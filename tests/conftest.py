from collections.abc import Iterator
from typing import TYPE_CHECKING

import pytest
from numpy.testing import assert_array_almost_equal
from rdkit import Chem
from rdkit.Chem.rdMolTransforms import ComputeCentroid

from prolif.datafiles import TOP, TRAJ, WATER_TOP, WATER_TRAJ, datapath
from prolif.interactions.base import _INTERACTIONS
from prolif.molecule import Molecule, sdf_supplier

# MDAnalysis support is removed - skip related functionality
pytest_plugins = []

try:
    from MDAnalysis import Universe
    from MDAnalysis.topology.guessers import guess_atom_element

    _HAS_MDANALYSIS = True
except ImportError:
    _HAS_MDANALYSIS = False
    Universe = None
    guess_atom_element = None

if TYPE_CHECKING:
    from prolif.molecule import BaseRDKitMol


def pytest_sessionstart(session: pytest.Session) -> None:
    if not datapath.exists():
        pytest.exit(
            f"Example data files are not accessible: {datapath!s} does not exist",
        )
    vina_path = datapath / "vina"
    if not vina_path.exists():
        pytest.exit(
            f"Example Vina data files are not accessible: {vina_path!s} does not exist",
        )
    # ugly patch to add Mixin class as attribute to pytest so that we don't have to
    # worry about relative imports in the test codebase
    pytest.BaseTestMixinRDKitMol = BaseTestMixinRDKitMol  # type: ignore[attr-defined]


# MDAnalysis-dependent fixtures are disabled
@pytest.fixture(scope="session")
def u():
    if not _HAS_MDANALYSIS:
        pytest.skip("MDAnalysis not available")
    return Universe(TOP, TRAJ)


@pytest.fixture(scope="session")
def rdkit_mol() -> Chem.Mol:
    return Chem.MolFromPDBFile(TOP, removeHs=False)


@pytest.fixture(scope="session")
def ligand_ag(u):
    if not _HAS_MDANALYSIS:
        pytest.skip("MDAnalysis not available")
    return u.select_atoms("resname LIG")


@pytest.fixture(scope="session")
def ligand_rdkit(ligand_ag) -> Chem.Mol:
    if not _HAS_MDANALYSIS:
        pytest.skip("MDAnalysis not available")
    return ligand_ag.convert_to.rdkit()  # type: ignore[no-any-return]


@pytest.fixture(scope="session")
def ligand_mol(ligand_ag) -> Molecule:
    if not _HAS_MDANALYSIS:
        pytest.skip("MDAnalysis not available")
    # Use SDF supplier instead of from_mda to avoid NotImplementedError
    from prolif.datafiles import datapath

    sdf_path = datapath / "vina" / "vina_output.sdf"
    from prolif.molecule import sdf_supplier

    mols = list(sdf_supplier(str(sdf_path)))
    return mols[0]  # Return first molecule from SDF


@pytest.fixture(scope="session")
def protein_ag(u, ligand_ag):
    if not _HAS_MDANALYSIS:
        pytest.skip("MDAnalysis not available")
    return u.select_atoms("protein and byres around 6.5 group ligand", ligand=ligand_ag)


@pytest.fixture(scope="session")
def protein_rdkit(protein_ag) -> Chem.Mol:
    if not _HAS_MDANALYSIS:
        pytest.skip("MDAnalysis not available")
    return protein_ag.convert_to.rdkit()  # type: ignore[no-any-return]


@pytest.fixture(scope="session")
def protein_mol(protein_ag) -> Molecule:
    if not _HAS_MDANALYSIS:
        pytest.skip("MDAnalysis not available")
    # Use PDB file instead of from_mda to avoid NotImplementedError
    from prolif.datafiles import datapath

    pdb_path = datapath / "top.pdb"
    from rdkit import Chem

    mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False)
    if mol is None:
        pytest.skip("Could not load protein PDB file")
    return Molecule.from_rdkit(mol)


@pytest.fixture(scope="session")
def sdf_suppl() -> sdf_supplier:
    path = str(datapath / "vina" / "vina_output.sdf")
    return sdf_supplier(path)


def from_mol2(filename: str) -> Molecule:
    # MOL2-based tests are deprecated due to MDAnalysis dependency removal
    # Skip all MOL2-based tests
    pytest.skip(
        f"MOL2 file {filename} tests skipped due to MDAnalysis dependency removal"
    )


@pytest.fixture(scope="session")
def benzene() -> Molecule:
    return from_mol2("benzene.mol2")


@pytest.fixture(scope="session")
def cation() -> Molecule:
    return from_mol2("cation.mol2")


@pytest.fixture(scope="session")
def cation_false() -> Molecule:
    return from_mol2("cation_false.mol2")


@pytest.fixture(scope="session")
def anion() -> Molecule:
    return from_mol2("anion.mol2")


@pytest.fixture(scope="session")
def ftf() -> Molecule:
    return from_mol2("facetoface.mol2")


@pytest.fixture(scope="session")
def etf() -> Molecule:
    return from_mol2("edgetoface.mol2")


@pytest.fixture(scope="session")
def chlorine() -> Molecule:
    return from_mol2("chlorine.mol2")


@pytest.fixture(scope="session")
def bromine() -> Molecule:
    return from_mol2("bromine.mol2")


@pytest.fixture(scope="session")
def hb_donor() -> Molecule:
    return from_mol2("donor.mol2")


@pytest.fixture(scope="session")
def hb_acceptor() -> Molecule:
    return from_mol2("acceptor.mol2")


@pytest.fixture(scope="session")
def hb_acceptor_false() -> Molecule:
    return from_mol2("acceptor_false.mol2")


@pytest.fixture(scope="session")
def xb_donor() -> Molecule:
    return from_mol2("xbond_donor.mol2")


@pytest.fixture(scope="session")
def xb_acceptor() -> Molecule:
    return from_mol2("xbond_acceptor.mol2")


@pytest.fixture(scope="session")
def xb_acceptor_false_xar() -> Molecule:
    return from_mol2("xbond_acceptor_false_xar.mol2")


@pytest.fixture(scope="session")
def xb_acceptor_false_axd() -> Molecule:
    return from_mol2("xbond_acceptor_false_axd.mol2")


@pytest.fixture(scope="session")
def ligand() -> Molecule:
    return from_mol2("ligand.mol2")


@pytest.fixture(scope="session")
def metal() -> Molecule:
    return from_mol2("metal.mol2")


@pytest.fixture(scope="session")
def metal_false() -> Molecule:
    return from_mol2("metal_false.mol2")


@pytest.fixture
def cleanup_dummy() -> Iterator[None]:
    yield
    _INTERACTIONS.pop("Dummy", None)


@pytest.fixture(scope="session")
def water_u():
    if not _HAS_MDANALYSIS:
        pytest.skip("MDAnalysis not available")
    return Universe(WATER_TOP, WATER_TRAJ)


@pytest.fixture(scope="session")
def water_atomgroups(water_u):
    if not _HAS_MDANALYSIS:
        pytest.skip("MDAnalysis not available")
    ligand = water_u.select_atoms("resname QNB")
    protein = water_u.select_atoms("protein and resid 399:404")
    water = water_u.select_atoms("segid WAT and (resid 17 or resid 83)")
    return ligand, protein, water


@pytest.fixture(scope="session")
def water_mols(water_atomgroups):
    if not _HAS_MDANALYSIS:
        pytest.skip("MDAnalysis not available")
    # Load from files instead of using from_mda to avoid NotImplementedError
    from prolif.datafiles import datapath

    # Load ligand from SDF
    sdf_path = datapath / "vina" / "vina_output.sdf"
    from prolif.molecule import sdf_supplier

    lig_mols = list(sdf_supplier(str(sdf_path)))
    lig_mol = lig_mols[0]

    # Load protein from PDB
    pdb_path = datapath / "water_m2.pdb"
    from rdkit import Chem

    prot_rdkit = Chem.MolFromPDBFile(str(pdb_path), removeHs=False)
    if prot_rdkit is None:
        pytest.skip("Could not load protein PDB for water tests")
    prot_mol = Molecule.from_rdkit(prot_rdkit)

    # Create a simple water molecule
    water_rdkit = Chem.MolFromSmiles("O")
    water_rdkit = Chem.AddHs(water_rdkit)
    from rdkit.Chem import rdDistGeom

    rdDistGeom.EmbedMolecule(water_rdkit)
    water_mol = Molecule.from_rdkit(water_rdkit)

    return lig_mol, prot_mol, water_mol


class BaseTestMixinRDKitMol:
    def test_init(self, mol: "BaseRDKitMol") -> None:
        assert isinstance(mol, Chem.Mol)

    def test_centroid(self, mol: "BaseRDKitMol") -> None:
        expected = ComputeCentroid(mol.GetConformer())
        assert_array_almost_equal(mol.centroid, expected)

    def test_xyz(self, mol: "BaseRDKitMol") -> None:
        expected = mol.GetConformer().GetPositions()
        assert_array_almost_equal(mol.xyz, expected)
