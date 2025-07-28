from collections.abc import Iterator
from typing import TYPE_CHECKING, cast

import matplotlib as mpl
import pytest
from matplotlib import pyplot as plt

import prolif as plf
from prolif.exceptions import RunRequiredError
from prolif.plotting.barcode import Barcode

try:
    import MDAnalysis as mda
    _HAS_MDANALYSIS = True
except ImportError:
    _HAS_MDANALYSIS = False
    pytest.skip("MDAnalysis not available", allow_module_level=True)

if TYPE_CHECKING:
    from MDAnalysis.core.groups import AtomGroup
    from MDAnalysis.core.universe import Universe


class TestBarcode:
    @pytest.fixture(scope="class", autouse=True)
    def backend(self) -> Iterator[None]:
        """Set backend to Agg to avoid TCL bugs with Windows."""
        backend = mpl.get_backend()
        mpl.use("Agg")
        yield
        mpl.use(backend)

    @pytest.fixture(scope="class")
    def simple_fp(self) -> plf.Fingerprint:
        return plf.Fingerprint(count=False)

    @pytest.fixture(scope="class")
    def count_fp(self) -> plf.Fingerprint:
        return plf.Fingerprint(count=True)

    @pytest.fixture(scope="class", params=["simple_fp", "count_fp"])
    def fp(self, request: pytest.FixtureRequest) -> plf.Fingerprint:
        return cast(plf.Fingerprint, request.getfixturevalue(request.param))

    @pytest.fixture(scope="class")
    def fp_run(
        self,
        ligand_mol: plf.Molecule,
        protein_mol: plf.Molecule,
        fp: plf.Fingerprint,
    ) -> plf.Fingerprint:
        # Use run_from_iterable to create proper IFP structure  
        fp.run_from_iterable([ligand_mol], protein_mol)
        return fp

    @pytest.fixture(scope="class")
    def barcode(self, fp_run: plf.Fingerprint) -> Barcode:
        return Barcode.from_fingerprint(fp_run)

    def test_display(self, barcode: Barcode) -> None:
        ax = barcode.display()
        assert isinstance(ax, plt.Axes)

    def test_display_kwargs(self, barcode: Barcode) -> None:
        ax = barcode.display(
            figsize=(1, 2),
            dpi=200,
            interactive=True,
            n_frame_ticks=2,
            residues_tick_location="bottom",
            xlabel="foobar",
            subplots_kwargs={},
            tight_layout_kwargs={"pad": 2},
        )
        assert isinstance(ax, plt.Axes)

    def test_fp_plot_barcode(self, fp_run: plf.Fingerprint) -> None:
        ax = fp_run.plot_barcode(
            figsize=(1, 2),
            dpi=200,
            interactive=True,
            n_frame_ticks=2,
            residues_tick_location="bottom",
            xlabel="foobar",
            subplots_kwargs={},
            tight_layout_kwargs={"pad": 2},
        )
        assert isinstance(ax, plt.Axes)

    def test_from_fingerprint_raises_not_executed(self) -> None:
        fp = plf.Fingerprint()
        with pytest.raises(
            RunRequiredError,
            match="Please run the fingerprint analysis before attempting to display"
            " results",
        ):
            Barcode.from_fingerprint(fp)
