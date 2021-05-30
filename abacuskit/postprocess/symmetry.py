'''
Date: 2021-05-12 20:44:55
LastEditors: jiyuyang
LastEditTime: 2021-05-12 20:44:55
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from collections import defaultdict
from typing import Counter, Dict, OrderedDict, Tuple, Union

import numpy as np
import seekpath
import spglib
from abacuskit.calculations.structure import Kpt, Stru
from abacuskit.utils.constants import BOHR_TO_A, CrySysNumber


class Spacegroup:
    """Interface to Spglib"""

    def __init__(self, stru: Stru = None):
        """Initialize Spacegroup object

        :params stru: Stru object
        """

        self.stru = stru

    @property
    def spgcell(self):
        """Set Spglib cell from Stru object"""

        lattice = np.array(self.stru.cell)*self.stru.lat0 * \
            BOHR_TO_A  # in Cartesian
        positions = []
        numbers = []
        magmoms = []
        for index, elem in enumerate(self.stru.elements):
            num = self.stru.numbers[elem]
            for n in range(num):
                numbers.append(index)
                # fractional atomic positions
                positions.append(self.stru.scaled_positions[elem][n])
                magmoms.append(self.stru.magmoms[elem])
        if sum(magmoms):
            spgcell = (lattice, positions, numbers, magmoms)
        else:
            spgcell = (lattice, positions, numbers)

        return spgcell

    def modified_stru(self, spgcell) -> Union[Stru, None]:
        """Modify Stru object based on spgcell. 

        :params spgcell: spglib cell (lattice, positions, numbers)
        """

        index = 0
        scaled_positions = defaultdict(list)
        numbers_dict = OrderedDict()
        if spgcell:
            lattice, positions, numbers = spgcell
        else:
            return None
        tot_numbers = Counter(numbers)
        for i, elem in enumerate(self.stru.pps):
            numbers_dict[elem] = tot_numbers[i]
            for j in range(tot_numbers[i]):
                scaled_positions[elem].append(positions[index+j])
            index += tot_numbers[i]
        lattice /= self.stru.lat0*BOHR_TO_A

        return Stru(self.stru.lat0, lattice, self.stru.pps, scaled_positions=scaled_positions, orbitals=self.stru.orbitals, masses=self.stru.masses, magmoms=self.stru.magmoms, move=self.stru.move, abfs=self.stru.abfs)

    def get_spacegroup(self, symprec: float = 1e-5, angle_tolerance: float = -1, symbol_type: int = 0) -> Union[str, None]:
        """Return International space group short symbol and number or Schoenflies symbol

        :params symprec: distance tolerance in Cartesian coordinates to find crystal symmetry. Default: 1e-5
        :params angle_tolerance: tolerance of angle between basis vectors in degrees to be tolerated in the symmetry finding. Default: -1
        :params symbol_type: 0 - International, 1 - Schoenflies
        """

        return spglib.get_spacegroup(self.spgcell, symprec, angle_tolerance, symbol_type)

    def get_symmetry(self, symprec: float = 1e-5, angle_tolerance: float = -1) -> Union[Dict[str, np.ndarray], None]:
        """Return dictionary of symmetry operations

        :params symprec: distance tolerance in Cartesian coordinates to find crystal symmetry. Default: 1e-5
        :params angle_tolerance: tolerance of angle between basis vectors in degrees to be tolerated in the symmetry finding. Default: -1
        :return dictionary with keys "rotation", "translation" and "equivalent_atoms"
        """

        return spglib.get_symmetry(self.spgcell, symprec, angle_tolerance)

    def refine_cell(self, symprec: float = 1e-5, angle_tolerance: float = -1) -> Union[Stru, None]:
        """Standardized crystal structure following space group type

        :params symprec: distance tolerance in Cartesian coordinates to find crystal symmetry. Default: 1e-5
        :params angle_tolerance: tolerance of angle between basis vectors in degrees to be tolerated in the symmetry finding. Default: -1
        :return lattice, atomic scaled positions, and atomic numbersthat are symmetrized following space group type
        """

        newcell = spglib.refine_cell(self.spgcell, symprec, angle_tolerance)

        return self.modified_stru(newcell)

    def find_primitive(self, symprec: float = 1e-5, angle_tolerance: float = -1) -> Union[Stru, None]:
        """Return a primitive cell

        :params symprec: distance tolerance in Cartesian coordinates to find crystal symmetry. Default: 1e-5
        :params angle_tolerance: tolerance of angle between basis vectors in degrees to be tolerated in the symmetry finding. Default: -1
        :return lattice, atomic scaled positions, and atomic numbersthat are symmetrized following space group type
        """

        newcell = spglib.find_primitive(self.spgcell, symprec, angle_tolerance)

        return self.modified_stru(newcell)

    def get_symmetry_dataset(self, symprec: float = 1e-5, angle_tolerance: float = -1, hall_number: int = 0) -> dict:
        """Search symmetry dataset from an input cell.

        :params symprec: distance tolerance in Cartesian coordinates to find crystal symmetry. Default: 1e-5
        :params angle_tolerance: tolerance of angle between basis vectors in degrees to be tolerated in the symmetry finding. Default: -1
        :params hall_number: If a serial number of Hall symbol from 1 to 530 is given, the database corresponding to the Hall symbol is made.
        """

        return spglib.get_symmetry_dataset(self.spgcell, symprec, angle_tolerance, hall_number)

    def get_kpath(self, numbers: list = [], with_time_reversal: bool = True, recipe: str = 'hpkot', threshold: float = 1e-7, symprec: float = 1e-5, angle_tolerance: float = -1) -> Tuple[Stru, Kpt]:
        """Obtain primitive cell and its band paths in the Brillouin zone of crystal structures.

        :params numbers: list of number of k points between two adjacent special k points.
        :params with_time_reversal: if time-reversal symmetry is present or not. Default: True
        :params recipe: currently only the string 'hpkot' is supported. Default: 'hpkot'
        :params threshold: numerical threshold. Default: 1e-7
        :params symprec: distance tolerance in Cartesian coordinates to find crystal symmetry. Default: 1e-5
        :params angle_tolerance: tolerance of angle between basis vectors in degrees to be tolerated in the symmetry finding. Default: -1
        """

        if numbers:
            set_already = True
        else:
            set_already = False

        res = seekpath.get_path(
            self.spgcell, with_time_reversal, recipe, threshold, symprec, angle_tolerance)
        newcell = (res['primitive_lattice'],
                   res['primitive_positions'], res['primitive_types'])
        stru = self.modified_stru(newcell)
        klabel = []
        special_k = []
        for k_tuple in res['path']:
            klabel.append(k_tuple[0])
            if not set_already:
                numbers.append(20)
        klabel.append(k_tuple[1])
        if not set_already:
            numbers.append(1)
        for label in klabel:
            special_k.append(res['point_coords'][label])
        kpt = Kpt(mode='Line', numbers=numbers,
                  special_k=special_k, klabel=klabel)

        return stru, kpt

    @staticmethod
    def get_crystal_system(number: int):
        """Return name of crystal system from space group number

        :params number: space group number
        """

        for key, value in CrySysNumber.items():
            if number in value:
                return key
