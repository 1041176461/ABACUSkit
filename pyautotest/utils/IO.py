'''
Date: 2021-05-16 16:43:23
LastEditors: jiyuyang
LastEditTime: 2021-05-16 16:43:24
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.utils.constants import BOHR_TO_A
from pyautotest.utils.typings import *
from pyautotest.utils.tools import get_input_line, list_elem2str, list_elem_2float, list_elem_2int, search_sentence, skip_notes
from pyautotest.calculations.structure import Stru, Kpt, Orb

import re
import json
from collections import defaultdict
from typing import Tuple

def read_json(filename: str_PathLike) -> dict:
    """ Read json file and return dict
    
    :params filename: json file
    """
    with open(filename, 'r') as file:
        text = json.load(file)
    return text

def write_json(filename: str_PathLike, new_filename: str_PathLike, **kwargs) -> str_PathLike:
    """ Read json file and modify some key-values, then write it to a new json file

    :params filename: json file to be read
    :params new_filename: json file to be written
    :params **kwargs: any key-value to be written to new_filename
    """
    with open(filename, 'r') as file:
        text = json.load(file)
        for key, value in kwargs.items():
            text[key] = value
    with open(new_filename, 'w') as file:
        json.dump(text, file, indent=4)
    return new_filename

def read_cif(filename: str_PathLike) -> Tuple[tuple, dict]:
    """Read cif file, return lattice and position
    
    :params filename: cif file
    """
    res = {}
    with open(filename, 'r') as file:
        for line in file:
            if re.search("_cell_length_a", line):
                a = float(line.split()[1])
            if re.search("_cell_length_b", line):
                b = float(line.split()[1])
            if re.search("_cell_length_c", line):
                c = float(line.split()[1])
            if re.search("_cell_angle_alpha", line):
                alpha = float(line.split()[1])
            if re.search("_cell_angle_beta", line):
                beta = float(line.split()[1])
            if re.search("_cell_angle_gamma", line):
                gamma = float(line.split()[1])
            if re.search("_atom_site_fract_z", line):
                position = defaultdict(list)
                for atom in file:
                    elem, x, y, z = atom.split()
                    position[elem].append((float(x), float(y), float(z)))
        lattice = (a, b, c, alpha, beta, gamma)
    
    return lattice, position

def write_input(input_dict:dict, filename:str_PathLike="INPUT"):
    """Write INPUT file based on input_dict
        
    :params input_dict: dict of input parameters
    """

    with open(filename, 'w') as file:
        file.write(get_input_line(input_dict))

def read_stru(ntype: int, stru_file: str_PathLike) -> Stru:
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
                cell.append(list_elem_2float(skip_notes(file.readline()).split()))

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
                    R_tmp.append(list_elem_2float(line.split()[:3]))
                    move_tmp.append(list_elem_2int(line.split()[3:]))
                if ctype == "Direct":
                    scaled_positions[elem] = R_tmp
                elif ctype == "Cartesian":
                    positions[elem] = R_tmp
                elif ctype == "Cartesian_angstrom":
                    positions_angstrom_lat0[elem] = R_tmp
                move[elem] = move_tmp

    return Stru(lat0, cell, pps, positions=positions, scaled_positions=scaled_positions, positions_angstrom_lat0=positions_angstrom_lat0, orbitals=orbitals, masses=masses, magmoms=magmoms, move=move, abfs=abfs)

def read_kpt(kpt_file: str_PathLike) -> Kpt:
    """Read k-points file
    
    :params kpt_file: absolute path of k-points file
    """
    number = 0
    with open(kpt_file,"r") as file:
        if search_sentence(file, ["K_POINTS", "KPOINTS", "K"]):
            number = int(file.readline())
        mode = search_sentence(file, ["Gamma", "MP", "Line"])
        if mode in ["Gamma", "MP"]:
            line = skip_notes(file.readline()).split()
            numbers = list_elem_2int(line[:3])
            offset = list_elem_2float(line[3:])
            return Kpt(mode, numbers, special_k=[], offset=offset)
        elif mode == "Line":
            special_k = []
            numbers = []
            klabel = []
            for k in range(number):
                line = file.readline()
                if re.match("#", line):
                    continue
                else:
                    linesplit = line.split(maxsplit=4)
                special_k.append(list_elem_2float(linesplit[:3]))
                numbers.append(int(linesplit[3]))
                if len(linesplit) == 5:
                    klabel.append(linesplit[4].strip('#\n '))
            return Kpt(mode, numbers, special_k, offset=[], klabel=klabel)

def read_orb(orbital: str_PathLike) -> Orb:
    """Read orbital file
    
    :params orbital: absolute path of orbital file
    """

    Nu = []
    with open(orbital, "r") as file:
        for line in file:
            line = skip_notes(line)
            if line.startswith("Element"):
                element = line.split()[-1]
            elif line.startswith("Energy Cutoff(Ry)"):
                ecut = float(line.split()[-1])
            elif line.startswith("Radius Cutoff(a.u.)"):
                rcut = float(line.split()[-1])
            elif line.startswith("Number of"):
                Nu.append(int(line.split()[-1]))
                
    return Orb(element=element, ecut=ecut, rcut=rcut, Nu=Nu, datafile=orbital)

def POSCAR_to_Stru(filename:str_PathLike="", pps:Dict_str_str={}, orbitals:Dict_str_str={}, masses:Dict_str_float={}, magmoms:Dict_str_float={}, move:Dict_str_int={}, abfs:Dict_str_str={}) -> Stru:
    """
    Convert POSCAR of VASP to STRU of ABACUS
    """

    def skip_notes(line):
        line = re.compile("[' ']+").sub(",", line).strip("\n").split(",")
        return line[1:]

    with open(filename, 'r') as f:
        tag = f.readline()
        lat0 = float(f.readline())/BOHR_TO_A
        cell = []
        for i in range(3):
            cell.append(list_elem_2float(skip_notes(f.readline())))
        elem = skip_notes(f.readline())
        num = skip_notes(f.readline())
        elem_num = dict(zip(elem, num))
        ctype = f.readline().strip()
        positions = {}
        scaled_positions = {}
        for j, k in elem_num.items():
            R_tmp = []
            for m in range(int(k)):
                R_tmp.append(list_elem_2float(skip_notes(f.readline())))
            if ctype == "Direct":
                scaled_positions[j] = R_tmp
            elif ctype == "Cartesian":
                positions[j] = R_tmp

    return Stru(lat0, cell, pps, positions=positions, scaled_positions=scaled_positions, orbitals=orbitals, masses=masses, magmoms=magmoms, move=move, abfs=abfs)