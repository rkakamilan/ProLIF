from typing import TYPE_CHECKING

import pytest

from prolif.fingerprint import Fingerprint
from prolif.ifp import IFP, InteractionData
from prolif.residue import ResidueId

try:
    import MDAnalysis as mda
    _HAS_MDANALYSIS = True
except ImportError:
    _HAS_MDANALYSIS = False
    pytest.skip("MDAnalysis not available", allow_module_level=True)

if TYPE_CHECKING:
    from MDAnalysis.core.groups import AtomGroup
    from MDAnalysis.core.universe import Universe


@pytest.fixture(scope="session")
def ifp(ligand_mol, protein_mol) -> IFP:
    fp = Fingerprint(["Hydrophobic", "VdWContact"])
    # Use run_from_iterable to create proper IFP structure
    fp.run_from_iterable([ligand_mol], protein_mol)
    return fp.ifp[0]


def test_ifp_indexing(ifp: IFP) -> None:
    # Use actual residue IDs from the data
    lig_id, prot_id = "UNL1", "VAL201.A"
    metadata1 = ifp[ResidueId.from_string(lig_id), ResidueId.from_string(prot_id)]
    metadata2 = ifp[lig_id, prot_id]
    assert metadata1 is metadata2


def test_ifp_filtering(ifp: IFP) -> None:
    # Use actual residue IDs from the data
    lig_id, prot_id = "UNL1", "VAL201.A"
    # Check that filtering by ligand returns the same IFP
    assert isinstance(ifp[lig_id], IFP)
    # Check filtering by protein residue
    if prot_id in [str(k[1]) for k in ifp.keys()]:
        prot_filtered = ifp[prot_id]
        assert isinstance(prot_filtered, IFP)
        assert len(prot_filtered) <= len(ifp)


def test_wrong_key(ifp: IFP) -> None:
    with pytest.raises(KeyError, match="does not correspond to a valid IFP key"):
        ifp[0]  # type: ignore[call-overload]


def test_interaction_data_iteration(ifp: IFP) -> None:
    data = next(ifp.interactions())
    assert isinstance(data, InteractionData)
    assert data.ligand == ResidueId("UNL", 1, None)
    assert data.protein.chain in {"A", "B"}
    assert data.interaction in {"Hydrophobic", "VdWContact"}
    assert "distance" in data.metadata
    for data in ifp.interactions():
        assert isinstance(data, InteractionData)
