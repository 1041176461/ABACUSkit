'''
Date: 2021-04-28 21:11:37
LastEditors: jiyuyang
LastEditTime: 2021-04-29 23:32:04
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.calculations.baseclass import JobCalculation
from pyautotest.utils.tools import folder_name, list_elem2str, read_json
from pyautotest.calculations.structure import Stru
from pyautotest.schedulers.data import Code
from pyautotest.utils.typings import *
from pyautotest.calculations.plugins.dis import distance

import os
import glob
import json
import textwrap
import shutil
import numpy as np
from pathlib import Path

class SetDimers(JobCalculation):
    """Set and calculate dimers for atomic orbital(LCAO)"""

    def __init__(self, element: str, ecutwfc: float, nbands: int, ref_band: int, Nu: list, rcut: float, pps: Dict_str_str, sigma: float=0.001, target: int=1, **kwargs) -> None:
        """Set input parameters of dimers calculation

        :params element: string of element name
        :params ecutwfc: energy cutoff for plane wave functions, the unit is Rydberg
        :params nbands: number of bands to calculate
        :params ref_band: number of reference bands
        :params Nu: size of LCAO, e.g. `[4,3,2,1]` means 4s3p2d1f
        :params rcut: cutoff of wavefunctions(a.u.)
        :params pps: dict of pseudopotential file
        :params sigma: energy range for smearing, the unit is Rydberg.
        :params target: fitting target, 0 - wave function, 1 - wave function and its gradient
        """

        super().__init__(**kwargs)
        self.element = element
        self.ref_band = ref_band
        self.Nu = Nu
        self.rcut = rcut
        self.pps = pps
        self.input_dict = {}
        self.input_dict["ntype"] = 1
        self.input_dict["nbands"] = nbands
        self.input_dict["ecutwfc"] = ecutwfc
        self.input_dict["basis_type"] = "pw"
        self.input_dict["wannier_card"] = "INPUTw"
        self.input_dict["calculation"] = "scf"
        self.input_dict["nspin"] = 1
        self.input_dict["lmaxmax"] = len(self.Nu)-1
        self.input_dict["symmetry"] = 0
        self.input_dict["dr2"] = 2.0e-8
        self.input_dict["niter"] = 1000
        self.input_dict["smearing"] = "gauss"
        self.input_dict["sigma"] = sigma
        self.input_dict["mixing_type"] = "pulay"
        self.input_dict["mixing_beta"] = 0.4
        self.input_dict["mixing_ndim"] = 8
        self.input_dict["printe"] = 1
        self.folder_list = []
        self.dimer_num = len(distance[self.element])
        self.target = target

    def _prepare(self, **kwargs):
        """Prepare input files for dimers calculation"""

        for i_dis in distance[self.element]:
            self.folder_list.append(self.create_input(i_dis))

    def _get_input_line(self):
        """Return input lines in INPUT file"""

        lines = []
        lines.append("INPUT_PARAMETERS")
        for key, value in self.input_dict.items():
            if value:
                lines.append(f"{key.ljust(30)}{value}")
        return '\n'.join(lines)

    def create_input(self, i_dis: float):
        """Create dimer directory and its input files

        :params i_dis: distance of dimers
        """

        folder = Path(folder_name(self.element, self.rcut, i_dis))
        folder.mkdir(parents=True, exist_ok=False)

        # Pseudopotential
        for i in self.pps.values():
            shutil.copy(i, folder)

        # INPUT
        input_lines = self._get_input_line()
        with open(folder/"INPUT", 'w') as file:
            file.write(input_lines)
        with open(folder/"INPUTw", "w") as file:
            file.write(textwrap.dedent(f"""\
                    WANNIER_PARAMETERS
                    rcut 10
                    out_spillage 2
                    """))
        with open(folder/"INPUTs", "w") as file:
            file.write(textwrap.dedent(f"""\
                    INPUT_ORBITAL_INFORMATION
                    <SPHERICAL_BESSEL>
                    0           // smooth or not
                    0.1         // sigma
                    {self.input_dict["ecutwfc"]}          // energy cutoff for spherical bessel functions(Ry)
                    {self.rcut}          // cutoff of wavefunctions(a.u.)
                    1.0e-12     // tolerence
                    </SPHERICAL_BESSEL>
                    """))
        
        # KPT
        with open(folder/"KPT", "w") as file:
            file.write(textwrap.dedent("""
                K_POINTS
                0
                Gamma
                1 1 1 0 0 0
            """))

        # STRU
        lat0=30
        cell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        pps = {self.element : self.pps[self.element]}
        masses = {self.element : 1}
        magmoms = {self.element : 0}
        if self.element in ["Na","Li","K","Ca"]:
            positions_angstrom_lat0 = {self.element : [[0.0, 0.0, 0.0], [0.0, 0.0, i_dis], [0.0, i_dis*0.86603, i_dis*0.5]]}
            move = {self.element : [[0, 0, 0], [0, 0, 0], [0, 0, 0]]}
        else:
            positions_angstrom_lat0 = {self.element : [[0.0, 0.0, 0.0], [0.0, 0.0, i_dis]]}
            move = {self.element : [[0, 0, 0], [0, 0, 0]]}
        obj = Stru(positions_angstrom_lat0=positions_angstrom_lat0, lat0=lat0, cell=cell, pps=pps, masses=masses, magmoms=magmoms, move=move)
        obj.write_stru(folder/"STRU")

        return str(folder)

    def _execute(self, command: Command, count: int=1, **kwargs):
        """Execute calculation

        :params command: string to execute a scf calculation
        :params count: how many times does it take for all dimers to be calculated, each time `(number of dimers)/count` be calculated
        """

        filename = "folders"
        origin = [os.path.join(folder, "test.0.dat") for folder in self.folder_list]
        if self.target == 1:
            linear = [[os.path.join(folder, "test.1.dat") for folder in self.folder_list]]
        data = {"origin" : origin, "linear" : linear}
        with open(filename, 'w') as file:  
            json.dump(data, file, indent=4)

        def set_commands(one_list, command):
            if isinstance(command, Code):
                command = command.run_line()
            commands = []
            for folder in one_list:
                commands.append(f"(cd {folder};{command};cd ../) &")
            line =  '\n'.join(commands)+"\nwait\n"
            return line

        current_path = Path.cwd()
        lines = []
        all_list = np.array_split(self.folder_list, count)
        for one_list in all_list:
            lines.append(set_commands(one_list, command))
        os.system('\n'.join(lines))
        os.chdir(current_path)

    def _check(self, **kwargs):
        """Check if job is finished"""
        
        for folder in self.folder_list:
            if not glob.glob(f"{folder}/test.*.dat"):
                raise FileNotFoundError(f"'test.*.dat' file not found in {folder}, may be something wrong.")

    def _parse(self, **kwargs) -> dict:
        """Parse output of a finished job and write a file named 'folder' for optimize ABFS orbitals
        
        :return: some result
        """

        res = {"dimer_matrix":1}

        return res

class OptLCAO(JobCalculation):
    """Optimize atomic orbitals(LCAO)"""

    def __init__(self, element: str, ecutwfc: float, ref_band: int, Nu: list, rcut: float, target: int=1, cal_T: bool=True, cal_smooth: bool=True, dr: float=0.01, lr: float=0.01, **kwargs):
        """Set input parameters of LCAOs optimization

        :params element: string of element name
        :params ecutwfc: energy cutoff for plane wave functions, the unit is Rydberg
        :params ref_band: number of reference bands
        :params Nu: size of LCAO, e.g. `[4,3,2,1]` means 4s3p2d1f
        :params target: fitting target, 0 - wave function, 1 - wave function and its gradient
        :params cal_T: if fit kinetic energy: Default: True
        :params cal_smooth: if use smooth method. Default: True
        """

        super().__init__(**kwargs)
        self.element = element
        self.ecutwfc = ecutwfc
        self.ref_band = ref_band
        self.Nu = Nu
        self.rcut = rcut
        self.dr = dr
        self.lr = lr
        self.target = target
        self.dimer_num = len(distance[self.element])
        self.cal_T = "true" if cal_T else "false"
        self.cal_smooth = "true" if cal_smooth else "false"

    def _prepare(self, **kwargs):
        """Prepare input files for optimizing ABFs e.g. input.json"""

        self.folder_opt = Path("opt_orb_"+"-".join(list_elem2str(self.Nu)))
        self.folder_opt.mkdir(exist_ok=False)

        if not Path("folders").exists():
            raise FileNotFoundError("'folders' which is a out file of `SetDimers` calculations not found.")

        with open(self.folder_opt/"input.json", 'w') as file:
            json.dump(self.set_input(), file, indent=4)

    def set_input(self):
        input_params = {
            "file_list": read_json("folders"),
            "info":
            {
                "Nt_all":       [self.element],
                "Nu":           {self.element:self.Nu},
                "Nb_true":	    [self.ref_band]*self.dimer_num,
                "weight":	    [1]*self.dimer_num,
                "Rcut":		    {self.element:self.rcut},
                "dr":		    {self.element:self.dr},
                "Ecut":		    {self.element:self.ecutwfc},
                "lr":		    self.lr,
                "cal_T":        self.cal_T,
                "cal_smooth":   self.cal_smooth
            },
            "C_init_info":
            { 
                "init_from_file": False 
            },
            "V_info":
            {
                "init_from_file":	True,
                "same_band":		True
            }
        }
        return input_params

    def _execute(self, command: Command, **kwargs):
        """Execute calculation"""
        
        if isinstance(command, Code):
            line = command.run_line()
        elif isinstance(command, str):
            line = command

        current_path = Path.cwd()
        os.chdir(self.folder_opt)
        os.environ["MKL_THREADING_LAYER"] = "GNU"
        os.system(line)
        os.chdir(current_path)

    def _check(self, **kwargs):
        """Check if job is finished"""

        filename = f"{self.folder_opt}/orb_{self.element}.dat"
        if not glob.glob(filename):
            raise FileNotFoundError(f"{filename} not found.")

    def _parse(self, **kwargs) -> dict:
        """Parse output of a finished job"""

        # Orbitals
        filename = f"{self.folder_opt}/orb_{self.element}.dat"
        if self.target == 0:
            datafile = f"{self.folder_opt}/psi_{self.element}.dat"
            name = "PSI"
        elif self.target == 1:
            datafile = f"{self.folder_opt}/dpsi_{self.element}.dat"
            name = "DPSI"
        os.rename(filename, datafile)

        res = {name:1}
        return res

class LCAO(JobCalculation):
    """Linear Combination of Atomic Orbitals"""

    def __init__(self, element: str, ecutwfc: float, nbands: int, ref_band: int, Nu: list, rcut: float, pps: Dict_str_str, sigma: float=0.001, target: int=1, cal_T: bool=True, cal_smooth: bool=True, dr: float=0.01, lr: float=0.01, **kwargs) -> None:
        """Set input parameters of dimers calculation

        :params element: string of element name
        :params ecutwfc: energy cutoff for plane wave functions, the unit is Rydberg
        :params nbands: number of bands to calculate
        :params ref_band: number of reference bands
        :params Nu: size of LCAO, e.g. `[4,3,2,1]` means 4s3p2d1f
        :params rcut: cutoff of wavefunctions(a.u.)
        :params pps: dict of pseudopotential file
        :params sigma: energy range for smearing, the unit is Rydberg.
        :params target: fitting target, 0 - wave function, 1 - wave function and its gradient
        :params cal_T: if fit kinetic energy: Default: True
        :params cal_smooth: if use smooth method. Default: True
        """

        super().__init__(**kwargs)
        self.element = element
        self.ecutwfc = ecutwfc
        self.nbands = nbands
        self.ref_band = ref_band
        self.Nu = Nu
        self.rcut = rcut
        self.pps = pps
        self.sigma = sigma
        self.target = target
        self.cal_T = cal_T
        self.cal_smooth = cal_smooth
        self.dr = dr
        self.lr = lr
    
    def _prepare(self, **kwargs):
        """Prepare input files for hybrid  e.g. input.json"""

        self.obj_setdimers = SetDimers(self.element, self.ecutwfc, self.nbands, self.ref_band, self.Nu, self.rcut, self.pps, self.sigma, self.target, **self.kwargs)
        self.obj_optlcao = OptLCAO(self.element, self.ecutwfc, self.ref_band, self.Nu, self.rcut, self.target, self.cal_T, self.cal_smooth, self.dr, self.lr, **self.kwargs)

    def _execute(self, command: Command, external_command: Code, **kwargs):
        """Execute calculation

        :params command: pyautotest.schedulers.data.Code` object to execute calculation
        """

        # set dimer
        self.obj_setdimers.calculate(command, **kwargs)

        # optimize LCAO
        self.obj_optlcao.calculate(external_command, **kwargs)

    def _check(self, **kwargs):
        """Check if job is finished"""
        pass

    def _parse(self, **kwargs) -> dict:
        """Parse output of a finished job
        
        :return: some result
        """

        res = {"LCAO":1}
        return res