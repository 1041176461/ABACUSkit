'''
Date: 2021-03-31 15:38:25
LastEditors: jiyuyang
LastEditTime: 2021-04-29 16:18:17
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.utils.tools import skip_notes, search_sentence, list_elem2str

import re
import numpy as np
import functools
import operator
import itertools
import decimal
from copy import deepcopy
from collections import defaultdict

def Direct2Cartesian(positions, cell):
    """Transform direct coordinates to Cartesian format in unit lat0
    
    :params postions: atomic direct coordinates in unit lat0
    :params cell: list of lattice vectors. If `format` is "Direct", it should be set. lattice vectors in ABACUS are scaled by the lattice constant. Default: []
    """

    for pos in positions:
        positions[pos] = np.dot(positions[pos], cell)

    return positions

def Cartesian_angstrom2Cartesian(positions):
    """Transform Cartesian_angstrom coordinates to Cartesian format in unit lat0
    
    :params postions: atomic Cartesian_angstrom coordinates in unit lat0
    """

    for pos in positions:
        positions[pos] = np.array(positions[pos])/0.529166

    return positions

def read_stru(ntype, stru_file=""):
    """
    Read `STRU` file

    :params stru_file: absolute path of `STRU` file
    """
    
    elements = []
    masses = {}
    pps = {}
    orbitals = {}
    cell = []
    magmoms = {}
    numbers = {}
    positions = defaultdict(list)
    scaled_positions = defaultdict(list)
    positions_angstrom_lat0 = defaultdict(list)
    move = defaultdict(list)
    abfs = {}
    with open(stru_file, "r") as file:
        if search_sentence(file, "ATOMIC_SPECIES"):
            for it in range(ntype):
                line = skip_notes(file.readline())
                elem, mass, pseudo = line.split()
                elements.append(elem)
                masses[elem] = float(mass)
                pps[elem] = pseudo

        if search_sentence(file, "NUMERICAL_ORBITAL"):
            for elem in elements:
                orbitals[elem] = skip_notes(file.readline())

        if search_sentence(file, "ABFS_ORBITALL"):
            for elem in elements:
                abfs[elem] = skip_notes(file.readline())

        if search_sentence(file, "LATTICE_CONSTANT"):
            lat0 = float(skip_notes(file.readline()).split()[0])

        if search_sentence(file, "LATTICE_VECTORS"):
            for i in range(3):
                cell.append(list(map(float, skip_notes(file.readline()).split())))

        if search_sentence(file, "ATOMIC_POSITIONS"):
            ctype = skip_notes(file.readline())

        for elem in elements:
            if search_sentence(file, elem):
                magmoms[elem] = float(skip_notes(file.readline()).split()[0])
                na = int(skip_notes(file.readline()).split()[0])
                numbers[elem] = na
                R_tmp = []
                move_tmp = []
                for i in range(na):
                    line = skip_notes(file.readline())
                    R_tmp.append(list(map(float, line.split()[:3])))
                    move_tmp.append(list(map(int, line.split()[3:])))
                if ctype == "Direct":
                    scaled_positions[elem] = R_tmp
                elif ctype == "Cartesian":
                    positions[elem] = R_tmp
                elif ctype == "Cartesian_angstrom":
                    positions_angstrom_lat0[elem] = R_tmp
                move[elem] = move_tmp

    return Stru(elements=elements, positions=positions, scaled_positions=scaled_positions, positions_angstrom_lat0=positions_angstrom_lat0, lat0=lat0, cell=cell, pps=pps, orbitals=orbitals, numbers=numbers, masses=masses, magmoms=magmoms, move=move, abfs=abfs)

class Stru:
    """ABACUS `STRU` file information"""

    def __init__(self, elements=[], positions={}, scaled_positions={}, positions_angstrom_lat0={}, lat0=None, cell=[], pps={}, orbitals={}, numbers={}, masses={}, magmoms={}, move={}, abfs={}) -> None:
        """Initialize Stru object
        
        :params elements: list, list of atomic species. Default: []
        :params positions: dict, key is element name and value is list of atomic Cartesian coordinates in unit lat0. Default: {}
        :params scaled_positions: dict, key is element name and value is list of atomic direct coordinates in unit lat0. This parameter had better not be set at the same time with `positions`. 
                                Because Cartesian coordinates will be calculated automatically based on direct coordinates and lattice vector `cell`. So self.positions still exists. Default: {}
        :params positions_angstrom_lat0: dict, key is element name and value is list of atomic Cartesian_angstrom coordinates in unit lat0. This parameter had better not be set at the same time with `positions`. 
                                Because Cartesian coordinates will be calculated automatically based on Cartesian_angstrom coordinates. So self.positions still exists. Default: {}
        :params lat0: int, lattice constant in unit bohr. Default: None
        :params cell: list, lattice vector in unit `lat0`. Default: []
        :params pps: dict, dict of pseudopotential file. Default: {}
        :params orbitals: dict, dict of orbital file. Default: {}
        :params numbers: dict, dict of numbers of atom. Default: {}
        :params masses: dict, dict of atomic mass. Default: {}
        :params magmoms: dict, dict of magnetic moment. Default: {}
        :params move: dict, key is element name and value is list of 1 or 0, `1` means atom can move. Default: {}
        :params abfs: dict, dict of ABFs for hybrid functional calculation. Default: {}
        """

        self.elements = elements
        self.ntype = len(self.elements)
        self.lat0 = lat0
        self.cell = cell
        if positions and scaled_positions and positions_angstrom_lat0:
            raise TypeError("'positions', 'scaled_positions' and `positions_angstrom_lat0` can not be set simultaneously")
        elif positions:
            self._ctype = "Cartesian"
            self.scaled_positions = scaled_positions
            self.positions = positions
        elif scaled_positions:
            self._ctype = "Direct"
            self.scaled_positions = scaled_positions
            self.positions = Direct2Cartesian(scaled_positions.copy(), self.cell)
        elif positions_angstrom_lat0:
            self._ctype = "Cartesian_angstrom"
            self.positions_angstrom_lat0 = positions_angstrom_lat0
            self.positions = Cartesian_angstrom2Cartesian(positions_angstrom_lat0)
        else:
            raise TypeError("One of 'positions' and 'scaled_positions' must be set")
        self.pps = pps
        self.orbitals = orbitals
        self.numbers = numbers
        if len(masses) == 0:
            self.masses = {elem:1 for elem in self.elements}
        else:
            self.masses = masses
        if len(magmoms) == 0:
            self.magmoms = {elem:0 for elem in self.elements}
        else:
            self.magmoms = magmoms
        if len(move) == 0:
            for i, elem in enumerate(self.elements):
                for j in range(self.numbers[i]):
                    move[elem].append([1, 1, 1])
        self.move = move
        self.abfs = abfs

    @property
    def positions_bohr(self):
        new_positions = deepcopy(self.positions)
        for pos in new_positions:
            new_positions[pos] = np.array(new_positions[pos]) * self.lat0

        return new_positions

    def get_stru(self):
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
                    line.append(" ".join(list_elem2str(self.positions[elem][j]))+" "+" ".join(list_elem2str(self.move[elem][j])))
                elif self._ctype == "Direct":
                    line.append(" ".join(list_elem2str(self.scaled_positions[elem][j]))+" "+" ".join(list_elem2str(self.move[elem][j])))
                elif self._ctype == "Cartesian_angstrom":
                    line.append(" ".join(list_elem2str(self.positions_angstrom_lat0[elem][j]))+" "+" ".join(list_elem2str(self.move[elem][j])))
            line.append(empty_line) 

        return '\n'.join(line)

    def write_stru(self, filename):
        """write `STRU` file
        
        :params filename: absolute path of `STRU` file
        """

        with open(filename, 'w') as file:
            file.write(self.get_stru())

    def supercell_positions(self, kpt):
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
                        R_new.append(R[pos]+np.dot(np.array([ix,iy,iz]), lat_vec))
            R[pos] = np.concatenate(R_new)

        return R

def cal_dis(positions={}):
    """Calculate distance between two atoms
    
    :params positions: dict, key is element name and value is list of atomic Cartesian coordinates in unit lat0
    :return dis[T1,T2] = {..., i_dis:num, ...}
    """

    dis = dict()
    for T1,T2 in itertools.combinations_with_replacement(positions, 2):
        dis_TT = defaultdict(int)
        for ia1,ia2 in itertools.product(positions[T1], positions[T2]):
            i_dis = np.linalg.norm(ia1-ia2)
            dis_TT[i_dis] += 1
        dis[T1,T2] = dict(dis_TT)

    return dis

def cut_dis(dis, Rcut):
    """Cut off distance between two atoms based on their orbitals cut-off radius
    
    :params dis: dict, distance between two atoms. Its format is `dis[T1,T2] = {..., i_dis:num, ...}`
    :params Rcut: dict, orbitals cut-off radius
    """

    for T1,T2 in dis:
        Rcut_sum = Rcut[T1]+Rcut[T2]
        dis[T1,T2] = { i_dis:num for i_dis,num in dis[T1,T2].items() if i_dis<Rcut_sum }

    return dis

def round_dis(dis, precision):
    """Handle the precision of distance value
    
    :params dis: dict, distance between two atoms. Its format is `dis[T1,T2] = {..., i_dis:num, ...}`
    :params precision: floating point precision
    """
    dis_round = dict()
    for T1,T2 in dis:
        dis_TT = defaultdict(int)
        for i_dis,num in dis[T1,T2].items():
            i_dis = float(decimal.Decimal(i_dis).quantize(decimal.Decimal(str(precision)), rounding=decimal.ROUND_HALF_UP))	
            dis_TT[i_dis] += num
        dis_round[T1,T2] = dict(dis_TT)

    return dis_round

def delete_zero(dis):
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

def read_kpt(kpt_file=""):
    """Read `KPT` file
    
    :params kpt_file: absolute path of `KPT` file
    """
    with open(kpt_file,"r") as file:
        search_sentence(file, ["Gamma", "MP"])
        line = skip_notes(file.readline()).split()
        return list(map(int, line[:3])), list(map(float, line[3:]))

def read_orb(orbital):
    """Read orbital file
    
    :params orbital: absolute path of orbital file
    """
    am = []
    with open(orbital,"r") as file:
        for line in file:
            line = skip_notes(line)
            if line.startswith("Element"):
                element = line.split()[-1]
            elif line.startswith("Energy Cutoff(Ry)"):
                ecut = float(line.split()[-1])
            elif line.startswith("Radius Cutoff(a.u.)"):
                rcut = float(line.split()[-1])
            elif line.startswith("Number of"):
                am.append(int(line.split()[-1]))
                
    return Orb(element=element, ecut=ecut, rcut=rcut, am=am, datafile=orbital)

class Orb:
    """Orbital information"""

    def __init__(self, element="", ecut=None, rcut=None, am=[], datafile="") -> None:
        """Initialize Orb object
        
        :params element: string of element name
        :params ecut: energy cutoff in unit Rydberg
        :params rcut: radius cutoff in unit bohr
        :params am: list of number of orbitals of each angular momentum
        :params datafile: absolute path of orbital file
        """

        self.element = element
        self.ecut = ecut
        self.rcut = rcut
        self.am = am
        self.datafile = datafile

    @property
    def total(self):
        """Return total number of orbitals"""
        return functools.reduce(operator.add, ((2*l+1)*n for l,n in enumerate(self.am)))

#TODO: add a function to plot orbital