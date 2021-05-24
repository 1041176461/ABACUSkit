'''
Date: 2021-05-12 20:44:55
LastEditors: jiyuyang
LastEditTime: 2021-05-12 20:44:55
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.calculations.structure import Stru
from pyautotest.utils.constants import BOHR_TO_A
from pyautotest.utils.typings import Dict_str_str, Dict_str_list, Dict_str_float, Dict_str_int

import spglib
import numpy as np
from collections import defaultdict
from typing import OrderedDict, Tuple, Union, Dict

class Spacegroup:
    """Interface to Spglib"""

    def __init__(self, stru:Stru=None):
        """Initialize Spacegroup object
        
        :params stru: Stru object
        """

        self.stru = stru
    
    @property
    def spgcell(self):
        """Set Spglib cell from Stru object"""

        lattice = np.array(self.stru.cell)*self.stru.lat0*BOHR_TO_A # in Cartesian
        positions = []
        numbers = []
        magmoms = []
        for elem in self.stru.elements:
            num = self.stru.numbers[elem]
            numbers.append(num)
            for n in range(num):
                positions.append(self.stru.scaled_positions[elem][n]) # fractional atomic positions
            magmoms.append(self.stru.magmoms[elem])
        if sum(magmoms):
            spgcell = (lattice, positions, numbers, magmoms)
        else:
            spgcell = (lattice, positions, numbers)
        
        return spgcell

    @staticmethod
    def modified_stru(spgcell, lat0: float, pps: Dict_str_str, scaled_positions: Dict_str_list={}, orbitals: Dict_str_str={}, masses: Dict_str_float={}, magmoms: Dict_str_float={}, move: Dict_str_int={}, abfs: Dict_str_str={}):
        """Modify Stru object based on spgcell. 
        
        :params spgcell: spglib cell (lattice, positions, numbers)
        """

        index = 0
        scaled_positions = defaultdict(list)
        numbers_dict = OrderedDict()
        lattice, positions, numbers = spgcell
        for i, elem in enumerate(pps):
            numbers_dict[elem] = numbers[i]
            for j in range(numbers[i]):
                scaled_positions[elem].append(positions[index+j])
            index += numbers[i]
        lattice /= lat0*BOHR_TO_A
        
        return Stru(lat0, lattice, pps, scaled_positions=scaled_positions, orbitals=orbitals, masses=masses, magmoms=magmoms, move=move, abfs=abfs)

    def get_spacegroup(self, symprec=1e-5, angle_tolerance=-1, symbol_type=0) -> Union[str, None]:
        """Return International space group short symbol and number or Schoenflies symbol
        
        :params symprec: distance tolerance in Cartesian coordinates to find crystal symmetry. Default: 1e-5
        :params angle_tolerance: tolerance of angle between basis vectors in degrees to be tolerated in the symmetry finding. Default: -1
        :params symbol_type: 0 - International, 1 - Schoenflies
        """
        
        return spglib.get_spacegroup(self.spgcell, symprec, angle_tolerance, symbol_type)

    def get_symmetry(self, symprec=1e-5, angle_tolerance=-1) -> Union[Dict[str, np.ndarray], None]:
        """Return dictionary of symmetry operations
        
        :params symprec: distance tolerance in Cartesian coordinates to find crystal symmetry. Default: 1e-5
        :params angle_tolerance: tolerance of angle between basis vectors in degrees to be tolerated in the symmetry finding. Default: -1
        :return dictionary with keys "rotation", "translation" and "equivalent_atoms"
        """

        return spglib.get_symmetry(self.spgcell, symprec, angle_tolerance)

    def refine_cell(self, symprec=1e-5, angle_tolerance=-1) -> Union[Tuple[np.ndarray, np.ndarray, np.ndarray], None]:
        """Standardized crystal structure following space group type
        
        :params symprec: distance tolerance in Cartesian coordinates to find crystal symmetry. Default: 1e-5
        :params angle_tolerance: tolerance of angle between basis vectors in degrees to be tolerated in the symmetry finding. Default: -1
        :return lattice, atomic scaled positions, and atomic numbersthat are symmetrized following space group type
        """

        return spglib.refine_cell(self.spgcell, symprec, angle_tolerance)

    def find_primitive(self, symprec=1e-5, angle_tolerance=-1) -> Union[Tuple[np.ndarray, np.ndarray, np.ndarray], None]:
        """Return a primitive cell
        
        :params symprec: distance tolerance in Cartesian coordinates to find crystal symmetry. Default: 1e-5
        :params angle_tolerance: tolerance of angle between basis vectors in degrees to be tolerated in the symmetry finding. Default: -1
        :return lattice, atomic scaled positions, and atomic numbersthat are symmetrized following space group type
        """

        return spglib.find_primitive(self.spgcell, symprec, angle_tolerance)

    def get_symmetry_dataset(self, symprec=1e-5, angle_tolerance=-1, hall_number=0) -> dict:
        """Search symmetry dataset from an input cell.

        :params symprec: distance tolerance in Cartesian coordinates to find crystal symmetry. Default: 1e-5
        :params angle_tolerance: tolerance of angle between basis vectors in degrees to be tolerated in the symmetry finding. Default: -1
        :params hall_number: If a serial number of Hall symbol from 1 to 530 is given, the database corresponding to the Hall symbol is made.
        """

        return spglib.get_symmetry_dataset(self.spgcell, symprec, angle_tolerance, hall_number)