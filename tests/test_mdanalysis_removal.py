"""
Tests for MDAnalysis removal - error handling and deprecation messages
"""

import pytest
from rdkit import Chem

from prolif.molecule import Molecule
from prolif.utils import select_over_trajectory


class TestMDAnalysisRemoval:
    """Test that MDAnalysis-dependent functionality raises appropriate errors."""

    def test_molecule_from_mda_raises_error(self):
        """Test that Molecule.from_mda() raises NotImplementedError."""
        with pytest.raises(NotImplementedError, match="MDAnalysis dependency"):
            Molecule.from_mda(None, "protein")

    def test_select_over_trajectory_raises_error(self):
        """Test that select_over_trajectory() raises NotImplementedError."""
        with pytest.raises(
            NotImplementedError, match="select_over_trajectory.*no longer supported"
        ):
            select_over_trajectory("protein", None)

    def test_pdbqt_supplier_raises_error(self):
        """Test that pdbqt_supplier raises error when MDAnalysis is not available."""
        from prolif.molecule import pdbqt_supplier

        template = Chem.MolFromSmiles("CCO")
        with pytest.raises(
            NotImplementedError, match="PDBQT file support requires MDAnalysis"
        ):
            # This should raise NotImplementedError immediately
            next(iter(pdbqt_supplier(["dummy.pdbqt"], template)))

    def test_vdw_radii_data_independent(self):
        """Test that VdW radii data is available without MDAnalysis."""
        from prolif.interactions.constants import get_vdw_radius

        # Test common elements
        assert get_vdw_radius("C") == 1.7
        assert get_vdw_radius("N") == 1.55
        assert get_vdw_radius("O") == 1.52
        assert get_vdw_radius("H") == 1.2

    def test_basic_molecule_functionality_works(self):
        """Test that basic molecule functionality works without MDAnalysis."""
        mol = Chem.MolFromSmiles("CCO")
        # Add hydrogens to get the expected count
        mol = Chem.AddHs(mol)
        prolif_mol = Molecule.from_rdkit(mol)

        assert prolif_mol.n_residues == 1
        assert len(prolif_mol.GetAtoms()) == 9  # Including hydrogens
        assert prolif_mol[0].resid.name == "UNL"
