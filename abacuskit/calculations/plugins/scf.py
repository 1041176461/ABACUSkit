'''
Date: 2021-03-08 09:20:29
LastEditors: jiyuyang
LastEditTime: 2021-05-17 09:37:31
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

import re
import typing

import numpy as np
from abacuskit.calculations.baseclass import ABACUSCalculation
from abacuskit.calculations.structure import Kpt, Stru
from abacuskit.utils.IO import read_stru


class SCF(ABACUSCalculation):
    """SCF calculation"""

    def __init__(self, input_dict: dict, stru: typing.Optional[Stru], kpt: typing.Optional[Kpt], **kwargs) -> None:
        """Set input parameters of scf calcultion

        :params input_dict: dict of input parameters
        :params stru: object of `abacuskit.calculations.structure.Stru`
        :params kpt: object of `abacuskit.calculations.structure.Kpt`
        """

        super().__init__(input_dict, stru, kpt, **kwargs)
        self.input_dict["calculation"] = "scf"
        self.logfile = "OUT.test/running_scf.log"

    def _parse(self, **kwargs) -> dict:
        """parse output of scf calculation"""

        res = {}
        self.stru.set_energy(self.logfile)
        res["etot"] = self.stru.energy
        with open(self.logfile, 'r') as file:
            for line in file:
                if re.search("The Ionic Phase", line):
                    res["ionic phase"] = float(line.split()[3])
                if re.search("Electronic Phase", line):
                    res["electronic phase"] = float(line.split()[2])

        return res


class RELAX(SCF):
    """Relax calculation"""

    def __init__(self, input_dict: dict, stru: typing.Optional[Stru], kpt: typing.Optional[Kpt], **kwargs) -> None:
        """Set input parameters of relax calcultion

        :params input_dict: dict of input parameters
        :params stru: object of `abacuskit.calculations.structure.Stru`
        :params kpt: object of `abacuskit.calculations.structure.Kpt`
        """

        super().__init__(input_dict, stru, kpt, **kwargs)
        self.input_dict["calculation"] = "relax"
        self.input_dict["force"] = "1"
        self.logfile = "OUT.test/running_relax.log"

    def _parse(self, **kwargs) -> dict:
        """parse output of scf calculation"""

        res = super()._parse(**kwargs)
        with open(self.logfile, 'r') as file:
            for line in file:
                if re.search(r"TOTAL ATOM NUMBER = [0-9]+", line):
                    natom = int(re.search("[0-9]+", line).group())
                    force = np.zeros((natom, 3))
                if re.search("TOTAL-FORCE \(eV/Angstrom\)", line):
                    for i in range(3):
                        file.readline()
                    for i in range(natom):
                        element, fx, fy, fz = file.readline().split()
                        force[i] = (float(fx), float(fy), float(fz))
                    res["force_mean"] = force.mean()

        obj = read_stru(self.input_dict["ntype"], "OUT.test/STRU_ION_D")
        plist = []
        for elem in obj.elements:
            plist = np.mean(obj.positions[elem])
        res["positions_mean"] = np.mean(plist)
        res["lattice_mean"] = np.mean(obj.cell)

        return res


class CELL_RELAX(RELAX):
    """Cell relax calculation"""

    def __init__(self, input_dict: dict, stru: typing.Optional[Stru], kpt: typing.Optional[Kpt], **kwargs) -> None:
        """Set input parameters of cell-relax calcultion

        :params input_dict: dict of input parameters
        :params stru: object of `abacuskit.calculations.structure.Stru`
        :params kpt: object of `abacuskit.calculations.structure.Kpt`
        """

        super().__init__(input_dict, stru, kpt, **kwargs)
        self.input_dict["calculation"] = "cell-relax"
        self.input_dict["stress"] = "1"
        self.logfile = "OUT.test/running_cell-relax.log"

    def _parse(self, **kwargs) -> dict:
        """parse output of cell-relax calculation"""

        res = super()._parse(**kwargs)
        with open(self.logfile, 'r') as file:
            for line in file:
                if re.search("TOTAL-STRESS \(KBAR\)", line):
                    stress = np.zeros((3, 3))
                    for i in range(3):
                        file.readline()
                    for i in range(3):
                        stress[i] = np.array(
                            file.readline().split(), dtype=np.float)
                    res["stress_mean"] = stress.mean()

        return res
