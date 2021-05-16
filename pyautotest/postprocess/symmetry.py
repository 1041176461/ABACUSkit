'''
Date: 2021-05-12 20:44:55
LastEditors: jiyuyang
LastEditTime: 2021-05-12 20:44:55
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from typing import Tuple, Union, Dict
from pyautotest.calculations.structure import Stru
from pyautotest.utils.constants import BOHR_TO_A

import spglib
import numpy as np

class Spacegroup:
    """Interface to Spglib"""

    def __init__(self, stru:Stru=None):
        """Set Spglib cell from Stru object
        
        :params stru: Stru object
        """

        lattice = stru.cell*stru.lat0*BOHR_TO_A # in Cartesian
        positions = []
        numbers = []
        magmoms = []
        for elem in stru.elements:
            positions.append(stru.scaled_positions[elem]) # fractional atomic positions
            numbers.append(stru.numbers[elem])
            magmoms.append(stru.magmoms[elem])
        if sum(magmoms):
            self.cell = (lattice, positions, numbers, magmoms)
        else:
            self.cell = (lattice, positions, numbers)

    def get_spacegroup(self, symprec=1e-5, angle_tolerance=-1, symbol_type=0) -> Union[str, None]:
        """Return International space group short symbol and number or Schoenflies symbol
        
        :params symprec: distance tolerance in Cartesian coordinates to find crystal symmetry. Default: 1e-5
        :params angle_tolerance: tolerance of angle between basis vectors in degrees to be tolerated in the symmetry finding. Default: -1
        :params symbol_type: 0 - International, 1 - Schoenflies
        """
        
        return spglib.get_spacegroup(self.cell, symprec, angle_tolerance, symbol_type)

    def get_symmetry(self, symprec=1e-5, angle_tolerance=-1) -> Union[Dict[str, np.ndarray], None]:
        """Return dictionary of symmetry operations
        
        :params symprec: distance tolerance in Cartesian coordinates to find crystal symmetry. Default: 1e-5
        :params angle_tolerance: tolerance of angle between basis vectors in degrees to be tolerated in the symmetry finding. Default: -1
        :return dictionary with keys "rotation", "translation" and "equivalent_atoms"
        """

        return spglib.get_symmetry(self.cell, symprec, angle_tolerance)

    def refine_cell(self, symprec=1e-5, angle_tolerance=-1) -> Union[Tuple[np.ndarray, np.ndarray, np.ndarray], None]:
        """Standardized crystal structure following space group type
        
        :params symprec: distance tolerance in Cartesian coordinates to find crystal symmetry. Default: 1e-5
        :params angle_tolerance: tolerance of angle between basis vectors in degrees to be tolerated in the symmetry finding. Default: -1
        :return lattice, atomic scaled positions, and atomic numbersthat are symmetrized following space group type
        """

        return spglib.refine_cell(self.cell, symprec, angle_tolerance)

    def find_primitive(self, symprec=1e-5, angle_tolerance=-1) -> Union[Tuple[np.ndarray, np.ndarray, np.ndarray], None]:
        """Return a primitive cell
        
        :params symprec: distance tolerance in Cartesian coordinates to find crystal symmetry. Default: 1e-5
        :params angle_tolerance: tolerance of angle between basis vectors in degrees to be tolerated in the symmetry finding. Default: -1
        :return lattice, atomic scaled positions, and atomic numbersthat are symmetrized following space group type
        """

        return spglib.find_primitive(self.cell, symprec, angle_tolerance)

    def get_symmetry_dataset(self, symprec=1e-5, angle_tolerance=-1, hall_number=0) -> dict:
        """Search symmetry dataset from an input cell.

        :params symprec: distance tolerance in Cartesian coordinates to find crystal symmetry. Default: 1e-5
        :params angle_tolerance: tolerance of angle between basis vectors in degrees to be tolerated in the symmetry finding. Default: -1
        :params hall_number: If a serial number of Hall symbol from 1 to 530 is given, the database corresponding to the Hall symbol is made.
        """

        return spglib.get_symmetry_dataset(self.cell, symprec, angle_tolerance, hall_number)