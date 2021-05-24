'''
Date: 2021-05-24 16:22:09
LastEditors: jiyuyang
LastEditTime: 2021-05-24 16:22:09
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.utils.IO import read_json
from pyautotest.calculations.structure import Stru
from pyautotest.utils.typings import *
from pyautotest.utils.tools import list_elem_2float
from pyautotest.utils.constants import BOHR_TO_A

import re
from pathlib import Path

class Convert:
    """Convert to ABACUS STRU file"""

    @classmethod
    def write(cls, filename:str_PathLike="", pps:Dict_str_str={}, orbitals:Dict_str_str={}, masses:Dict_str_float={}, magmoms:Dict_str_float={}, move:Dict_str_int={}, abfs:Dict_str_str={}) -> None:
        """"Write STRU file based on `filename` or its suffix"""

        if Path(filename).name == "POSCAR" or Path(filename).suffix in (".vasp", ".poscar", ".VASP"):
            stru = cls.POSCAR_to_STRU(filename, pps, orbitals, masses, magmoms, move, abfs)
        else:
            raise NotImplementedError("This format are not supported now.")
        stru.write_stru("STRU")

    @classmethod
    def POSCAR_to_STRU(cls, filename:str_PathLike="", pps:Dict_str_str={}, orbitals:Dict_str_str={}, masses:Dict_str_float={}, magmoms:Dict_str_float={}, move:Dict_str_int={}, abfs:Dict_str_str={}) -> Stru:
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

    @classmethod
    def convert_cmdline(cls, args):
        if args.stru:
            text = read_json(args.stru)
            filename = text["filename"]
            pps = text["pps"]
            orbitals = text.pop("orbitals", {})
            masses = text.pop("masses", {})
            magmoms = text.pop("magmoms", {})
            move = text.pop("move", {})
            abfs = text.pop("abfs", {})
            Convert().write(filename, pps, orbitals, masses, magmoms, move, abfs)