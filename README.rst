ProLIF
======

.. list-table::
    :widths: 12 35

    * - **Documentation**
      - |docs|
    * - **CI**
      - |tests| |codecov| |codeql|
    * - **Builds**
      - |conda-version| |pypi-version| |build|
    * - **Dependencies**
      - |rdkit|
    * - **License**
      - |license|

Description
-----------

ProLIF (*Protein-Ligand Interaction Fingerprints*) is a tool designed to generate
interaction fingerprints for complexes made of ligands, protein, DNA or RNA molecules
extracted from docking simulations and experimental structures.

.. note::
   
   **MDAnalysis Dependency Removed**: As of this version, MDAnalysis dependency has been 
   removed to resolve LGPL3 licensing issues. The package now focuses on docking results 
   analysis and no longer supports MD trajectory analysis. This enables use under more 
   permissive licenses.

Supported Features
------------------

✅ **Currently Supported**:

* Protein-ligand interaction analysis from docking results
* PDB protein structure support  
* SDF ligand structure support
* MOL2 file support (limited via RDKit)
* All core interaction types (hydrogen bonds, hydrophobic, π-stacking, etc.)
* Interaction fingerprint generation and analysis
* Export to DataFrames and bitvectors

❌ **No Longer Supported**:

* MD trajectory analysis (XTC, TRR, DCD formats)
* PDBQT file support (requires conversion to PDB/SDF)
* ``Molecule.from_mda()`` method
* ``select_over_trajectory()`` function  
* Advanced water bridge analysis requiring trajectory data

Migration Guide
---------------

**Old Code**::

    import MDAnalysis as mda
    import prolif

    u = mda.Universe("protein.pdb", "trajectory.xtc") 
    ligand = prolif.Molecule.from_mda(u.select_atoms("resname LIG"))

**New Code**::

    from rdkit import Chem
    import prolif

    protein = Chem.MolFromPDBFile("protein.pdb", removeHs=False)
    protein_mol = prolif.Molecule.from_rdkit(protein)
    ligands = prolif.sdf_supplier("ligands.sdf")

Documentation
-------------

The installation instructions, documentation and tutorials can be found online on
`ReadTheDocs <https://prolif.readthedocs.io>`_.

Issues
------

If you have found a bug, please open an issue on the
`GitHub Issues <https://github.com/chemosim-lab/ProLIF/issues>`_ page.

Discussion
----------

If you have questions on how to use ProLIF, or if you want to give feedback or share
ideas and new features, please head to the
`GitHub Discussions <https://github.com/chemosim-lab/ProLIF/discussions>`_ page.

Contributing
------------

If you are interested in contributing to the project, please read the
`Contribution Guidelines <CONTRIBUTING.md>`_.

Citing ProLIF
-------------

Please refer to the `citation page <https://prolif.readthedocs.io/en/latest/source/citation.html>`_
on the documentation.

License
-------

Unless otherwise noted, all files in this directory and all subdirectories are
distributed under the Apache License, Version 2.0 ::

    Copyright 2017-2025 Cédric BOUYSSET

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.


.. |pypi-version| image:: https://img.shields.io/pypi/v/prolif.svg
   :target: https://pypi.python.org/pypi/prolif
   :alt: Pypi Version

.. |conda-version| image:: https://img.shields.io/conda/vn/conda-forge/prolif.svg
    :target: https://anaconda.org/conda-forge/prolif
    :alt: Conda-forge version

.. |build| image:: https://github.com/chemosim-lab/ProLIF/workflows/build/badge.svg
    :target: https://github.com/chemosim-lab/ProLIF/actions/workflows/build.yml
    :alt: Build status

.. |tests| image:: https://github.com/chemosim-lab/ProLIF/workflows/tests/badge.svg?branch=master
    :target: https://github.com/chemosim-lab/ProLIF/actions/workflows/ci.yml
    :alt: Tests status

.. |codecov| image:: https://codecov.io/gh/chemosim-lab/ProLIF/branch/master/graph/badge.svg?token=2FCHV08G8A
    :target: https://codecov.io/gh/chemosim-lab/ProLIF

.. |docs| image:: https://readthedocs.org/projects/prolif/badge/?version=latest
    :target: https://prolif.readthedocs.io/
    :alt: Documentation Status

.. |codeql| image:: https://github.com/chemosim-lab/ProLIF/workflows/CodeQL/badge.svg?branch=master
    :target: https://github.com/chemosim-lab/ProLIF/actions/workflows/codeql.yml
    :alt: Code quality

.. |license| image:: https://img.shields.io/pypi/l/prolif
    :target: http://www.apache.org/licenses/LICENSE-2.0
    :alt: License


.. |rdkit| image:: https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC
      :alt: Powered by RDKit
      :target: https://www.rdkit.org/
