'''
Date: 2021-03-31 15:38:25
LastEditors: jiyuyang
LastEditTime: 2021-04-29 16:18:17
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

import decimal
import functools
import itertools
import math
import operator
import re
import typing
from collections import OrderedDict, defaultdict
from copy import deepcopy
from io import TextIOWrapper

import numpy as np
from abacuskit.utils.constants import BOHR_TO_A
from abacuskit.utils.tools import list_elem2str
from abacuskit.utils.typings import *


def conventional_cell(la: float, lb: float, lc: float, alpha: float, beta: float, gamma: float, crystal_system: str):
    """Set up Bravais lattice vectors of the conventional cell

    :params la, lb, lc, alpha, beta, gamma: lattice parameters used to generate the cell
    :params crystal_system: crystal system: 'cubic', 'hexagonal', 'tetragonal', 'trigonal', 'triclinic', 'monoclinic', 'unknown'
    """

    coa = lc/la
    boa = lb/la
    alphar = alpha*math.pi/180
    betar = beta*math.pi/180
    gammar = gamma*math.pi/180
    if crystal_system == 'cubic':
        cell = np.identity(3)
    elif crystal_system == 'hexagonal':
        cell = np.array([[math.sin(gammar), math.cos(gammar), 0],
                         [0, 1, 0],
                         [0, 0, coa]])
    elif crystal_system == 'tetragonal' or crystal_system == 'orthorhombic':
        cell = np.array([[1, 0, 0],
                         [0, boa, 0],
                         [0, 0, coa]])
    elif crystal_system == 'trigonal':
        # Hexagonal cell taken as conventional
        if not abs(gamma-120) < 0.0002:
            gammar = 120*math.pi/180
        cell = np.array([[math.sin(gammar), math.cos(gammar), 0],
                         [0, 1, 0],
                         [0, 0, coa]])
    elif crystal_system == 'triclinic' or crystal_system == 'monoclinic' or crystal_system == 'unknown':
        angfac1 = (math.cos(alphar) - math.cos(betar)
                   * math.cos(gammar))/math.sin(gammar)
        angfac2 = math.sqrt(math.sin(gammar)**2 - math.cos(betar)**2 - math.cos(alphar)**2
                            + 2*math.cos(alphar)*math.cos(betar)*math.cos(gammar))/math.sin(gammar)
        cell = np.array([[1, 0, 0],
                         [boa*math.cos(gammar), boa *
                          math.sin(gammar), 0],
                         [coa*math.cos(betar), coa*angfac1, coa*angfac2]])
    else:
        raise ValueError("No support for "+crystal_system+" crystal systems.")

    return cell*la

# STRU


def Direct2Cartesian(positions: Dict_str_list, cell: Dict_str_float) -> Dict_str_list:
    """Transform direct coordinates to Cartesian format in unit lat0

    :params postions: atomic direct coordinates in unit lat0
    :params cell: list of lattice vectors.
    """

    new_positions = deepcopy(positions)
    for pos in new_positions:
        new_positions[pos] = np.dot(new_positions[pos], cell)

    return new_positions


def Cartesian2Direct(positions: Dict_str_list, cell: Dict_str_float) -> Dict_str_list:
    """Transform Cartesian coordinates to direct format in unit lat0

    :params postions: atomic Cartesian coordinates in unit lat0
    :params cell: list of lattice vectors.
    """

    new_positions = deepcopy(positions)
    for pos in new_positions:
        new_positions[pos] = np.dot(new_positions[pos], np.linalg.inv(cell))

    return new_positions


def Cartesian_angstrom2Cartesian(positions: Dict_str_list) -> Dict_str_list:
    """Transform Cartesian_angstrom coordinates to Cartesian format in unit lat0

    :params postions: atomic Cartesian_angstrom coordinates in unit lat0
    """

    new_positions = deepcopy(positions)
    for pos in new_positions:
        new_positions[pos] = np.array(new_positions[pos])/BOHR_TO_A

    return new_positions


class Stru:
    """ABACUS `STRU` file information"""

    def __init__(self, lat0: float=1/BOHR_TO_A, cell: list=[], pps: Dict_str_str={}, positions: Dict_str_list = {}, scaled_positions: Dict_str_list = {}, positions_angstrom_lat0: Dict_str_list = {}, orbitals: Dict_str_str = {}, masses: Dict_str_float = {}, magmoms: Dict_str_float = {}, move: Dict_str_int = {}, abfs: Dict_str_str = {}) -> None:
        """Initialize Stru object

        :params lat0: float, lattice constant in unit bohr.
        :params cell: list, lattice vector in unit `lat0`.
        :params pps: dict, dict of pseudopotential file.
        :params positions: dict, key is element name and value is list of atomic Cartesian coordinates in unit lat0.
        :params scaled_positions: dict, key is element name and value is list of atomic direct coordinates in unit lat0. This parameter had better not be set at the same time with `positions`.
                                Because Cartesian coordinates will be calculated automatically based on direct coordinates and lattice vector `cell`. So self.positions still exists.
        :params positions_angstrom_lat0: dict, key is element name and value is list of atomic Cartesian_angstrom coordinates in unit lat0. This parameter had better not be set at the same time with `positions`.
                                Because Cartesian coordinates will be calculated automatically based on Cartesian_angstrom coordinates. So self.positions still exists.
        :params orbitals: dict, dict of orbital file.
        :params masses: dict, dict of atomic mass.
        :params magmoms: dict, dict of magnetic moment.
        :params move: dict, key is element name and value is list of 1 or 0, `1` means atom can move.
        :params abfs: dict, dict of ABFs for hybrid functional calculation.
        """

        self.lat0 = lat0
        self.cell = cell
        if positions and scaled_positions and positions_angstrom_lat0:
            raise TypeError(
                "'positions', 'scaled_positions' and `positions_angstrom_lat0` can not be set simultaneously")
        elif positions:
            self._ctype = "Cartesian"
            self.positions = positions
            self.scaled_positions = Cartesian2Direct(positions, self.cell)
        elif scaled_positions:
            self._ctype = "Direct"
            self.scaled_positions = scaled_positions
            self.positions = Direct2Cartesian(scaled_positions, self.cell)
        elif positions_angstrom_lat0:
            self._ctype = "Cartesian_angstrom"
            self.positions_angstrom_lat0 = positions_angstrom_lat0
            self.positions = Cartesian_angstrom2Cartesian(
                positions_angstrom_lat0)
        else:
            raise TypeError(
                "One of 'positions' and 'scaled_positions' must be set")
        self.elements = []
        self.numbers = OrderedDict()
        for elem in self.positions:
            self.elements.append(elem)
            self.numbers[elem] = len(self.positions[elem])
        self.pps = pps
        self.orbitals = orbitals
        if masses:
            self.masses = masses
        else:
            self.masses = {elem: 1 for elem in self.elements}
        if magmoms:
            self.magmoms = magmoms
        else:
            self.magmoms = {elem: 0 for elem in self.elements}
        if not move:
            move = defaultdict(list)
            for elem in self.elements:
                for j in range(self.numbers[elem]):
                    move[elem].append([1, 1, 1])
        self.move = move
        self.abfs = abfs

        self.energy = None
        self.efermi = None

    @property
    def positions_bohr(self):
        new_positions = deepcopy(self.positions)
        for pos in new_positions:
            new_positions[pos] = np.array(new_positions[pos]) * self.lat0

        return new_positions

    def get_stru(self) -> str:
        """Return the `STRU` file as a string"""

        empty_line = ''
        line = []
        line.append("ATOMIC_SPECIES")
        for elem in self.elements:
            line.append(f"{elem}\t{self.masses[elem]}\t{self.pps[elem]}")
        line.append(empty_line)

        if self.orbitals:
            line.append("NUMERICAL_ORBITAL")
            for elem in self.elements:
                line.append(f"{self.orbitals[elem]}")
            line.append(empty_line)

        if self.abfs:
            line.append("ABFS_ORBITAL")
            for elem in self.elements:
                line.append(f"{self.abfs[elem]}")
            line.append(empty_line)

        line.append("LATTICE_CONSTANT")
        line.append(str(self.lat0))
        line.append(empty_line)

        line.append("LATTICE_VECTORS")
        for i in range(3):
            line.append(" ".join(list_elem2str(self.cell[i])))
        line.append(empty_line)

        line.append("ATOMIC_POSITIONS")
        line.append(self._ctype)
        line.append(empty_line)
        for elem in self.elements:
            line.append(f"{elem}\n{self.magmoms[elem]}\n{self.numbers[elem]}")
            for j in range(self.numbers[elem]):
                if self._ctype == "Cartesian":
                    line.append(" ".join(list_elem2str(
                        self.positions[elem][j]))+" "+" ".join(list_elem2str(self.move[elem][j])))
                elif self._ctype == "Direct":
                    line.append(" ".join(list_elem2str(
                        self.scaled_positions[elem][j]))+" "+" ".join(list_elem2str(self.move[elem][j])))
                elif self._ctype == "Cartesian_angstrom":
                    line.append(" ".join(list_elem2str(
                        self.positions_angstrom_lat0[elem][j]))+" "+" ".join(list_elem2str(self.move[elem][j])))
            line.append(empty_line)

        return '\n'.join(line)

    def write_stru(self, filename: str = "STRU") -> None:
        """write `STRU` file

        :params filename: absolute path of `STRU` file
        """

        with open(filename, 'w') as file:
            file.write(self.get_stru())

    def get_cif(self, find_symmetry: bool = False):
        """ Get CIF format string

        :params find_symmetry: if use Spglib to find symmetry. Default: False
        """

        lattice = np.array(self.cell)*self.lat0 * BOHR_TO_A  # in Cartesian
        positions = []
        numbers = []
        magmoms = []
        label = []
        type_symbol = []
        for index, elem in enumerate(self.elements):
            num = self.numbers[elem]
            for n in range(num):
                numbers.append(index)
                label.append(f"{elem}{n+1}")
                type_symbol.append(f"{elem}")
                # fractional atomic positions
                positions.append(self.scaled_positions[elem][n])
                magmoms.append(self.magmoms[elem])
        if sum(magmoms):
            spgcell = (lattice, positions, numbers, magmoms)
        else:
            spgcell = (lattice, positions, numbers)

        from CifFile import CifBlock, CifFile
        cf = CifFile()
        cb = CifBlock()
        la = np.linalg.norm(lattice[0])
        lb = np.linalg.norm(lattice[1])
        lc = np.linalg.norm(lattice[2])
        cb['_cell_length_a'] = la
        cb['_cell_length_b'] = lb
        cb['_cell_length_c'] = lc
        cb['_cell_angle_alpha'] = math.acos(
            np.dot(self.cell[1], self.cell[2])/(lb*lc))/math.pi*180
        cb['_cell_angle_beta'] = math.acos(
            np.dot(self.cell[0], self.cell[2])/(la*lc))/math.pi*180
        cb['_cell_angle_gamma'] = math.acos(
            np.dot(self.cell[0], self.cell[1])/(la*lb))/math.pi*180

        from abacuskit.postprocess.symmetry import Spacegroup
        if find_symmetry:
            from spglib import get_symmetry_dataset
            dataset = get_symmetry_dataset(spgcell)
            if dataset:
                cb['_atom_site_label'] = [self.elements[i] for i in dataset['std_types']]
                cb['_atom_site_type_symbol'] = [f"{self.elements[i]}{i}" for i in dataset['std_types']]
                number = dataset['number']
                cb['_symmetry_cell_setting'] = Spacegroup(
                ).get_crystal_system(number)
                cb['_space_group_IT_number'] = number
                cb['_space_group_name_H-M_alt'] = dataset['international']
                cb['_space_group_name_Hall'] = dataset['hall']
                x, y, z = np.split(dataset["std_positions"], 3, axis=1)
            else:
                cb['_atom_site_label'] = label
                cb['_atom_site_type_symbol'] = type_symbol
                cb['_symmetry_cell_setting'] = Spacegroup().get_crystal_system(1)
                cb['_space_group_IT_number'] = 1
                cb['_space_group_name_H-M_alt'] = "P1"
                cb['_space_group_name_Hall'] = "P 1"
                x, y, z = np.split(np.array(positions), 3, axis=1)
        else:
            cb['_atom_site_label'] = label
            cb['_atom_site_type_symbol'] = type_symbol
            cb['_symmetry_cell_setting'] = Spacegroup().get_crystal_system(1)
            cb['_space_group_IT_number'] = 1
            cb['_space_group_name_H-M_alt'] = "P1"
            cb['_space_group_name_Hall'] = "P 1"
            x, y, z = np.split(np.array(positions), 3, axis=1)

        cb['_atom_site_fract_x'], cb['_atom_site_fract_y'], cb['_atom_site_fract_z'] = x.flatten(
        ), y.flatten(), z.flatten()
        cb.CreateLoop(['_atom_site_label', '_atom_site_type_symbol',
                       '_atom_site_fract_x', '_atom_site_fract_y', '_atom_site_fract_z'])

        cf['stru_to_cif'] = cb

        return str(cf)

    def write_cif(self, filename: str = "STRU.cif", fing_symmetry: bool = False):
        """write CIF file

        :params filename: absolute path of CIF file
        :params find_symmetry: if use Spglib to find symmetry. Default: False
        """

        with open(filename, 'w') as file:
            file.write(self.get_cif(fing_symmetry))

    def supercell_positions(self, kpt: list) -> Dict_str_list:
        """Return supercell atomic positions

        :params kpt: list of number of k-points in each directions
        """

        lat_vec = np.array(self.cell) * self.lat0
        R = deepcopy(self.positions_bohr)

        nx, ny, nz = kpt
        for pos in R:
            R_new = []
            for ix in range(nx):
                for iy in range(ny):
                    for iz in range(nz):
                        R_new.append(
                            R[pos]+np.dot(np.array([ix, iy, iz]), lat_vec))
            R[pos] = np.concatenate(R_new)

        return R

    @staticmethod
    def positions_dict2list(positions:dict):
        """Convert dict of positions to list
        
        :params positions: dict of positions, its keys are elements name and its values are positions
        """        
        new_positions = []
        for elem in positions:
            for pos in positions[elem]:
                new_positions.append(pos)

        return np.array(new_positions, dtype=float)

    @property
    def volume_bohr(self):
        return np.linalg.det(self.cell)*pow(self.lat0, 3)

    @property
    def volume(self):
        return self.volume_bohr*pow(BOHR_TO_A, 3)

    @classmethod
    def read_energy_from_file(cls, file: typing.Union[TextIOWrapper, str_PathLike], name: str):
        """Read energy from ABACUS running log file"""

        def read(f):
            for line in f:
                if re.search(name, line):
                    energy = float(line.split()[2])
            return energy

        if isinstance(file, TextIOWrapper):
            energy = read(file)
        elif isinstance(file, (str, PathLike)):
            with open(file, 'r') as f:
                energy = read(f)

        return energy

    def set_energy(self, file: typing.Union[TextIOWrapper, str_PathLike]):
        """Set energy"""

        self.energy = self.read_energy_from_file(file, "E_KohnSham")

    def set_efermi(self, file: TextIOWrapper):
        """Set fermi level"""

        self.efermi = self.read_energy_from_file(file, "E_Fermi")

    # TODO: add set_force, set_stress


def cal_dis(positions: dict, supercell_positions: dict) -> Dict_Tuple_Dict:
    """Calculate distance between two atoms

    :params positions: dict, key is element name and value is list of atomic Cartesian coordinates in unit lat0
    :params supercell_positions: dict, key is element name and value is supercell atomic positions
    :return dis[T1,T2] = {..., i_dis:num, ...}
    """

    dis = dict()
    for T1, T2 in itertools.combinations_with_replacement(positions, 2):
        dis_TT = defaultdict(int)
        for ia1, ia2 in itertools.product(positions[T1], supercell_positions[T2]):
            i_dis = np.linalg.norm(ia1-ia2)
            dis_TT[i_dis] += 1
        dis[T1, T2] = dict(dis_TT)

    return dis


def cut_dis(dis: Dict_Tuple_Dict, Rcut: Dict_str_float) -> Dict_Tuple_Dict:
    """Cut off distance between two atoms based on their orbitals cut-off radius

    :params dis: dict, distance between two atoms. Its format is `dis[T1,T2] = {..., i_dis:num, ...}`
    :params Rcut: dict, orbitals cut-off radius
    """

    for T1, T2 in dis:
        Rcut_sum = Rcut[T1]+Rcut[T2]
        dis[T1, T2] = {i_dis: num for i_dis,
                       num in dis[T1, T2].items() if i_dis < Rcut_sum}

    return dis


def round_dis(dis: Dict_Tuple_Dict, precision: float) -> Dict_Tuple_Dict:
    """Handle the precision of distance value

    :params dis: dict, distance between two atoms. Its format is `dis[T1,T2] = {..., i_dis:num, ...}`
    :params precision: floating point precision
    """
    dis_round = dict()
    for T1, T2 in dis:
        dis_TT = defaultdict(int)
        for i_dis, num in dis[T1, T2].items():
            i_dis = float(decimal.Decimal(i_dis).quantize(
                decimal.Decimal(str(precision)), rounding=decimal.ROUND_HALF_UP))
            dis_TT[i_dis] += num
        dis_round[T1, T2] = dict(dis_TT)

    return dis_round


def delete_zero(dis: Dict_Tuple_Dict) -> Dict_Tuple_Dict:
    """Delete zero in dict, list or set"""

    dis_new = dis.copy()
    while 0 in dis_new:
        if isinstance(dis_new, (list, set)):
            dis_new.remove(0)
        elif isinstance(dis_new, dict):
            dis_new.pop(0)
        else:
            raise TypeError

    return dis_new


# KPT
class Kpt:
    """K-points information"""

    def __init__(self, mode: str, numbers: list = [], special_k: list = [], offset: list = [0.0, 0.0, 0.0], klabel: list = []) -> None:
        """Initialize Kpt object

        :params mode: ‘Gamma’, ‘MP’ or 'Line'
        :params numbers: for ‘Gamma’ and ‘MP’ mode, list of three integers, for 'Line' mode, list of number of k points between two adjacent special k points.
        :params special_k: list of special k-points
        :params offset: offset of the k grid. Default: [0.0, 0.0, 0.0]
        """

        self.mode = mode
        self.numbers = numbers
        self.special_k = special_k
        self.offset = offset
        self.klabel = klabel

    def get_kpt(self) -> str:
        """Return the `KPT` file as a string"""

        line = []
        line.append("K_POINTS")
        if self.mode in ["Gamma", "MP"]:
            line.append("0")
            line.append(self.mode)
            line.append(" ".join(list_elem2str(self.numbers+self.offset)))
        elif self.mode == "Line":
            line.append(str(len(self.special_k)))
            line.append(self.mode)
            for i, k in enumerate(self.special_k):
                if self.klabel:
                    line.append(" ".join(list_elem2str(
                        k+[self.numbers[i]]))+'\t#'+self.klabel[i])
                else:
                    line.append(" ".join(list_elem2str(k+[self.numbers[i]])))

        return '\n'.join(line)

    def write_kpt(self, filename: str):
        """Write k-points file

        :params filename: absolute path of k-points file
        """

        with open(filename, 'w') as file:
            file.write(self.get_kpt())

    @property
    def full_kpath(self):
        """Uniform k-point path"""

        total_k = np.sum(self.numbers)
        spec_k_coor = np.array(self.special_k)
        interval = (np.roll(spec_k_coor, -1, axis=0) - spec_k_coor) / \
            np.reshape(self.numbers, (-1, 1))
        max_num = np.max(self.numbers)
        len_num = len(self.numbers)
        k_coor_span = np.zeros((len_num, max_num), dtype=np.float)
        X, Y, Z = np.split(spec_k_coor, 3, axis=1)
        i_X, i_Y, i_Z = np.split(interval, 3, axis=1)
        for i, j in enumerate(self.numbers):
            k_coor_span[i][:j] = np.arange(j)
        X = (i_X * k_coor_span + X.repeat(max_num, axis=1)).flatten()
        Y = (i_Y * k_coor_span + Y.repeat(max_num, axis=1)).flatten()
        Z = (i_Z * k_coor_span + Z.repeat(max_num, axis=1)).flatten()
        k_direct_coor = np.empty((3, total_k), dtype=np.float)
        k_direct_coor[0] = X[:total_k]
        k_direct_coor[1] = Y[:total_k]
        k_direct_coor[2] = Z[:total_k]

        return k_direct_coor.T

    @property
    def label_special_k(self):
        """Label special k-points based on `numbers` list"""

        index = np.cumsum(np.concatenate(([1], self.numbers), axis=0))[
            :len(self.special_k)]
        if self.klabel:
            return zip(self.klabel, index)
        else:
            return index

# Orb


class Orb:
    """Orbital information"""

    def __init__(self, element: str, ecut: float, rcut: float, Nu: list, datafile: str_PathLike) -> None:
        """Initialize Orb object

        :params element: string of element name
        :params ecut: energy cutoff in unit Rydberg
        :params rcut: radius cutoff in unit bohr
        :params Nu: list of number of orbitals of each angular momentum
        :params datafile: absolute path of orbital file
        """

        self.element = element
        self.ecut = ecut
        self.rcut = rcut
        self.Nu = Nu
        self.datafile = datafile

    @property
    def total(self):
        """Return total number of orbitals"""
        return functools.reduce(operator.add, ((2*l+1)*n for l, n in enumerate(self.Nu)))
# TODO: add a function to plot orbital
