'''
Date: 2021-05-16 16:43:23
LastEditors: jiyuyang
LastEditTime: 2021-05-16 16:43:24
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

import json
import re
from collections import defaultdict
from typing import List

import numpy as np
from abacuskit.calculations.structure import Kpt, Orb, Stru, conventional_cell
from abacuskit.postprocess.symmetry import Spacegroup
from abacuskit.utils.constants import BOHR_TO_A, Hall2Number, HM2Hall
from abacuskit.utils.tools import (get_input_line,
                                   list_elem_2float, list_elem_2int,
                                   search_sentence, skip_notes)
from abacuskit.utils.typings import *


def read_cif(filename: str_PathLike, pps: muti_Dict = [{}], orbitals: muti_Dict = [{}], masses: muti_Dict = [{}], magmoms: muti_Dict = [{}], move: muti_Dict = [{}], abfs: muti_Dict = [{}]) -> List[Stru]:
    """Read CIF(Crystallographic Information Framework) file and return `Stru` object

    :params filename: name of cif file
    """

    from CifFile import ReadCif
    strulist = []

    def list_elem2split(l):
        def list_split(s):
            return s.split('(')[0]
        return list(map(list_split, l))

    def get_elements(string):
        return re.search(r'([A-Z][a-z]?)', string).group(0)

    cf = ReadCif(filename)
    for index, name in enumerate(cf.keys()):
        data = cf[name]
        fract_x = list_elem2split(data["_atom_site_fract_x"])
        fract_y = list_elem2split(data["_atom_site_fract_y"])
        fract_z = list_elem2split(data["_atom_site_fract_z"])
        if "_atom_site_type_symbol" in data.keys():
            elements = list(map(get_elements, data["_atom_site_type_symbol"]))
        else:
            elements = list(map(get_elements, data["_atom_site_label"]))
        scaled_positions = defaultdict(list)
        for i, elem in enumerate(elements):
            scaled_positions[elem].append([np.array(fract_x[i], dtype=float), np.array(
                fract_y[i], dtype=float), np.array(fract_z[i], dtype=float)])

        la = float(data["_cell_length_a"].split('(')[0])
        lb = float(data["_cell_length_b"].split('(')[0])
        lc = float(data["_cell_length_c"].split('(')[0])
        alpha = float(data["_cell_angle_alpha"])
        beta = float(data["_cell_angle_beta"])
        gamma = float(data["_cell_angle_gamma"])

        if "_symmetry_cell_setting" in data.keys():
            crystal_system = data["_symmetry_cell_setting"]
        elif "_symmetry_Int_Tables_number" in data.keys():
            crystal_system = Spacegroup().get_crystal_system(
                data["_symmetry_Int_Tables_number"])
        elif "_space_group_name_Hall" in data.keys():
            crystal_system = Spacegroup().get_crystal_system(
                Hall2Number[data["_space_group_name_Hall"]])
        elif "_symmetry_space_group_name_H-M" in data.keys():
            crystal_system = Spacegroup().get_crystal_system(
                Hall2Number[HM2Hall[data["_symmetry_space_group_name_H-M "]]])
        else:
            crystal_system = 'unknown'

        cell = conventional_cell(
            la, lb, lc, alpha, beta, gamma, crystal_system)

        strulist.append(Stru(lat0=1/BOHR_TO_A, cell=cell, pps=pps[index], scaled_positions=scaled_positions,
                             orbitals=orbitals[index], masses=masses[index], magmoms=magmoms[index], move=move[index], abfs=abfs[index]))

    return strulist


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


def write_input(input_dict: dict, filename: str_PathLike = "INPUT"):
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
                cell.append(list_elem_2float(
                    skip_notes(file.readline()).split()))

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
    with open(kpt_file, "r") as file:
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
