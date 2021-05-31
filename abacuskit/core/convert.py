'''
Date: 2021-05-24 16:22:09
LastEditors: jiyuyang
LastEditTime: 2021-05-24 16:22:09
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

import re

from abacuskit.calculations.structure import Stru
from abacuskit.utils.constants import BOHR_TO_A
from abacuskit.utils.IO import read_cif, read_json, read_stru
from abacuskit.utils.tools import list_elem2str, list_elem_2float
from abacuskit.utils.typings import *


class Convert:
    """Convert file from one format to another"""

    @classmethod
    def POSCAR_to_STRU(cls, srcfile: str_PathLike, tofile: str_PathLike = "STRU", pps: Dict_str_str = {}, orbitals: Dict_str_str = {}, masses: Dict_str_float = {}, magmoms: Dict_str_float = {}, move: Dict_str_int = {}, abfs: Dict_str_str = {}) -> None:
        """Convert POSCAR of VASP to STRU of ABACUS

        :params srcfile: str_PathLike, source file
        :params tofile: str_PathLike, converted file
        :params pps: dict, dict of pseudopotential file. 
        :params orbitals: dict, dict of orbital file. 
        :params masses: dict, dict of atomic mass. 
        :params magmoms: dict, dict of magnetic moment. 
        :params move: dict, key is element name and value is list of 1 or 0, `1` means atom can move. 
        :params abfs: dict, dict of ABFs for hybrid functional calculation.
        """

        def skip_notes(line):
            line = re.compile("[' ']+").sub(",", line).strip("\n").split(",")
            return line[1:]

        with open(srcfile, 'r') as f:
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

        stru = Stru(lat0, cell, pps, positions=positions, scaled_positions=scaled_positions,
                    orbitals=orbitals, masses=masses, magmoms=magmoms, move=move, abfs=abfs)
        stru.write_stru(tofile)

    @classmethod
    def STRU_to_POSCAR(cls, srcfile: str_PathLike, tofile: str_PathLike, ntype: int) -> None:
        """Convert STRU of ABACUS to POSCAR of VASP

        :params srcfile: str_PathLike, source file
        :params tofile: str_PathLike, converted file
        """

        def convert_pos(elements: list, positions: dict):
            newpos = []
            for elem in elements:
                for pos in positions[elem]:
                    newpos.append(' '+''.ljust(20).join(list_elem2str(pos)))
            return newpos

        stru = read_stru(ntype, srcfile)
        line = [f"From {srcfile}"]
        line.append(str(stru.lat0*BOHR_TO_A))
        for i in stru.cell:
            line.append('\t'+'\t'.join(list_elem2str(i)))
        line.append(' '+''.ljust(10).join(stru.elements))
        line.append(' '+''.ljust(10).join(list_elem2str(stru.numbers.values())))
        line.append(stru._ctype)
        if stru._ctype == "Direct":
            line.append('\n'.join(convert_pos(
                stru.elements, stru.scaled_positions)))
        if stru._ctype == "Cartesian":
            line.append('\n'.join(convert_pos(stru.elements, stru.positions)))
        else:
            ValueError(f"Type {stru._ctype} is not supported for VASP.")
        with open(tofile, 'w') as file:
            file.write('\n'.join(line))

    @classmethod
    def CIF_to_STRU(cls, srcfile: str_PathLike, tofile: str_PathLike = "STRU", pps: Dict_str_str = {}, orbitals: Dict_str_str = {}, masses: Dict_str_float = {}, magmoms: Dict_str_float = {}, move: Dict_str_int = {}, abfs: Dict_str_str = {}) -> None:
        """Convert CIF(Crystallographic Information Framework) to STRU of ABACUS

        :params srcfile: str_PathLike, source file
        :params tofile: str_PathLike, converted file
        :params pps: dict, dict of pseudopotential file. 
        :params orbitals: dict, dict of orbital file. 
        :params masses: dict, dict of atomic mass. 
        :params magmoms: dict, dict of magnetic moment. 
        :params move: dict, key is element name and value is list of 1 or 0, `1` means atom can move. 
        :params abfs: dict, dict of ABFs for hybrid functional calculation.
        """

        stru = read_cif(srcfile, [pps], [orbitals], [
                        masses], [magmoms], [move], [abfs])[0]
        stru.write_stru(tofile)

    @classmethod
    def STRU_to_CIF(cls, srcfile: str_PathLike, tofile: str_PathLike, ntype: int, find_symmetry: bool = False) -> None:
        """Convert STRU of ABACUS to CIF(Crystallographic Information Framework)

        :params srcfile: str_PathLike, source file
        :params tofile: str_PathLike, converted file
        """

        stru = read_stru(ntype, srcfile)
        stru.write_cif(tofile, find_symmetry)

    @classmethod
    def _init_for_convert_STRU(cls, text):
        pps = text.pop("pps", {})
        orbitals = text.pop("orbitals", {})
        masses = text.pop("masses", {})
        magmoms = text.pop("magmoms", {})
        move = text.pop("move", {})
        abfs = text.pop("abfs", {})
        return pps, orbitals, masses, magmoms, move, abfs

    @classmethod
    def convert_cmdline(cls, args):
        if args.file:
            text = read_json(args.file)
            srcfile = text["srcfile"]
            convert = text["convert"]
            if convert == "poscar_to_stru":
                pps, orbitals, masses, magmoms, move, abfs = cls._init_for_convert_STRU(
                    text)
                tofile = text.pop("tofile", "STRU")
                Convert().POSCAR_to_STRU(srcfile, tofile, pps,
                                         orbitals, masses, magmoms, move, abfs)
            elif convert == "cif_to_stru":
                pps, orbitals, masses, magmoms, move, abfs = cls._init_for_convert_STRU(
                    text)
                tofile = text.pop("tofile", "STRU")
                Convert().CIF_to_STRU(srcfile, tofile, pps,
                                      orbitals, masses, magmoms, move, abfs)
            elif convert == "stru_to_poscar":
                ntype = text["ntype"]
                tofile = text.pop("tofile", "POSCAR")
                Convert().STRU_to_POSCAR(srcfile, tofile, ntype)
            elif convert == "stru_to_cif":
                ntype = text["ntype"]
                tofile = text.pop("tofile", "STRU.cif")
                Convert().STRU_to_CIF(srcfile, tofile, ntype)
