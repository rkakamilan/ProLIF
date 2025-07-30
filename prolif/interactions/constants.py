"""
Constants used by the package for interactions.
"""

from rdkit.Chem import GetPeriodicTable

# MDAnalysis-compatible VdW radii data (originally from MDAnalysis.topology.tables)
# Source: Standard atomic radii from various literature sources
_MDANALYSIS_VDWRADII = {
    "H": 1.1,
    "HE": 1.4,
    "LI": 1.82,
    "BE": 1.53,
    "B": 1.92,
    "C": 1.7,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "NE": 1.54,
    "NA": 2.27,
    "MG": 1.73,
    "AL": 1.84,
    "SI": 2.1,
    "P": 1.8,
    "S": 1.8,
    "CL": 1.75,
    "AR": 1.88,
    "K": 2.75,
    "CA": 2.31,
    "NI": 1.63,
    "CU": 1.4,
    "ZN": 1.39,
    "GA": 1.87,
    "GE": 2.11,
    "AA": 1.85,
    "SE": 1.9,
    "BR": 1.85,
    "KR": 2.02,
    "RR": 3.03,
    "SR": 2.49,
    "PD": 1.63,
    "AG": 1.72,
    "CD": 1.58,
    "IN": 1.93,
    "SN": 2.17,
    "SB": 2.06,
    "TE": 2.06,
    "I": 1.98,
    "XE": 2.16,
    "CS": 3.43,
    "BA": 2.68,
    "PT": 1.75,
    "AU": 1.66,
    "HH": 1.55,
    "TL": 1.96,
    "PB": 2.02,
    "BI": 2.07,
    "PO": 1.97,
    "AT": 2.02,
    "RN": 2.2,
    "FR": 3.48,
    "RA": 2.83,
    "U": 1.86,
}

VDWRADII: dict[str, float] = {
    symbol.capitalize(): radius for symbol, radius in _MDANALYSIS_VDWRADII.items()
}


def get_vdw_radius(element: str) -> float:
    """Get van der Waals radius for an element.

    Parameters
    ----------
    element : str
        Element symbol (e.g., 'C', 'N', 'O')

    Returns
    -------
    float
        van der Waals radius in Angstroms
    """
    return VDWRADII.get(element, 2.0)  # Default radius if not found


PT = GetPeriodicTable()
RDKIT_VDWRADII: dict[str, float] = {
    PT.GetElementSymbol(i): PT.GetRvdw(i) for i in range(1, 119)
}

# Table 1 of https://doi.org/10.1039/C3DT50599E
CSD_VDWRADII: dict[str, float] = {
    "H": 1.2,
    "He": 1.43,
    "Li": 2.12,
    "Be": 1.98,
    "B": 1.91,
    "C": 1.77,
    "N": 1.66,
    "O": 1.5,
    "F": 1.46,
    "Ne": 1.58,
    "Na": 2.5,
    "Mg": 2.51,
    "Al": 2.25,
    "Si": 2.19,
    "P": 1.9,
    "S": 1.89,
    "Cl": 1.82,
    "Ar": 1.83,
    "K": 2.73,
    "Ca": 2.62,
    "Sc": 2.58,
    "Ti": 2.46,
    "V": 2.42,
    "Cr": 2.45,
    "Mn": 2.45,
    "Fe": 2.44,
    "Co": 2.4,
    "Ni": 2.4,
    "Cu": 2.38,
    "Zn": 2.39,
    "Ga": 2.32,
    "Ge": 2.29,
    "As": 1.88,
    "Se": 1.82,
    "Br": 1.86,
    "Kr": 2.25,
    "Rb": 3.21,
    "Sr": 2.84,
    "Y": 2.75,
    "Zr": 2.52,
    "Nb": 2.56,
    "Mo": 2.45,
    "Tc": 2.44,
    "Ru": 2.46,
    "Rh": 2.44,
    "Pd": 2.15,
    "Ag": 2.53,
    "Cd": 2.49,
    "In": 2.43,
    "Sn": 2.42,
    "Sb": 2.47,
    "Te": 1.99,
    "I": 2.04,
    "Xe": 2.06,
    "Cs": 3.48,
    "Ba": 3.03,
    "La": 2.98,
    "Ce": 2.88,
    "Pr": 2.92,
    "Nd": 2.95,
    "Sm": 2.9,
    "Eu": 2.87,
    "Gd": 2.83,
    "Tb": 2.79,
    "Dy": 2.87,
    "Ho": 2.81,
    "Er": 2.83,
    "Tm": 2.79,
    "Yb": 2.8,
    "Lu": 2.74,
    "Hf": 2.63,
    "Ta": 2.53,
    "W": 2.57,
    "Re": 2.49,
    "Os": 2.48,
    "Ir": 2.41,
    "Pt": 2.29,
    "Au": 2.32,
    "Hg": 2.45,
    "Tl": 2.47,
    "Pb": 2.6,
    "Bi": 2.54,
    "Ac": 2.8,
    "Th": 2.93,
    "Pa": 2.88,
    "U": 2.71,  # noqa
    "Np": 2.82,
    "Pu": 2.81,
    "Am": 2.83,
    "Cm": 3.05,
    "Bk": 3.4,
    "Cf": 3.05,
    "Es": 2.7,
}

VDW_PRESETS = {
    "mdanalysis": VDWRADII,
    "rdkit": RDKIT_VDWRADII,
    "csd": CSD_VDWRADII,
}
