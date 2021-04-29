'''
Date: 2021-04-28 21:11:37
LastEditors: jiyuyang
LastEditTime: 2021-04-29 16:30:01
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.calculations.baseclass import ABACUSCalculation
from pyautotest.utils.tools import folder_name, list_elem2str, read_json
from pyautotest.calculations.structure import Stru
from pyautotest.schedulers.data import Code
from pyautotest.calculations.plugins.dis import distance

import os
import glob
import json
import textwrap
import shutil
import numpy as np
from pathlib import Path

class SetDimers(ABACUSCalculation):
    """Set and calculate dimers for atomic orbital(LCAO)"""

    def __init__(self, input_dict, src, element, Nu, rcut, ref_band, target=1, **kwargs) -> None:
        """Set input parameters of dimers calculation

        :params input_dict: dict of input parameters
        :params src: path of example which will be tested
        :params element: string of element name
        :params Nu: size of LCAO, e.g. `[4,3,2,1]` means 4s3p2d1f
        :params rcut: cutoff of wavefunctions(a.u.)
        :params ref_band: number of  reference bands
        :params target: fitting target, 0 - wave function, 1 - wave function and its gradient
        """

        super().__init__(input_dict, src, **kwargs)
        if "nbands" not in self.input_dict.keys():
            raise KeyError("`nbands` must be set in `input_dict`.")
        if "ecutwfc" not in self.input_dict.keys():
            raise KeyError("`ecutwfc` must be set in `input_dict`.")
        if "sigma" not in self.input_dict.keys():
            raise KeyError("`sigma` must be set in `input_dict`.")
        if not self.pps:
            raise TypeError("`pps` must be set.")
        self.element = element
        self.Nu = Nu
        self.input_dict["basis_type"] = "pw"
        self.input_dict["wannier_card"] = "INPUTw"
        self.input_dict["calculation"] = "scf"
        self.input_dict["nspin"] = 1
        self.input_dict["lmaxmax"] = len(self.Nu)-1
        self.input_dict["symmetry"] = 0
        self.input_dict["dr2"] = 2.0e-8
        self.input_dict["niter"] = 1000
        self.input_dict["smearing"] = "gauss"
        self.input_dict["mixing_type"] = "pulay"
        self.input_dict["mixing_beta"] = 0.4
        self.input_dict["mixing_ndim"] = 8
        self.input_dict["printe"] = 1
        self.rcut = rcut
        self.ref_band = ref_band
        self.folder_list = []
        self.dimer_num = len(distance[self.element])
        self.target = target

    def _prepare(self, dst, **kwargs):
        """Prepare input files for dimers calculation e.g. INPUT, STRU, KPT, orbital and pseudopotential files
        
        :params dst: path of working directory
        """

        for i_dis in distance[self.element]:
            self.folder_list.append(self.create_input(i_dis, dst))

    def create_input(self, i_dis, dst):
        """Create dimer directory and its input files
        
        :params T1: type 1
        :params T2: type 2
        :params i_dis: distance between `T1` and `T2`
        :params dst: directory where to create subdirectory named `T1`-`T2`_`i_dis`
        """

        folder = Path(dst, folder_name(self.element, self.rcut, i_dis))
        folder.mkdir(parents=True, exist_ok=False)

        # Pseudopotential
        for i in self.pps.values():
            shutil.copyfile(Path(self.src, i), Path(folder, i))

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
                    1           // smooth or not
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
            numbers = {self.element : 3}
            positions_angstrom_lat0 = {self.element : [[0.0, 0.0, 0.0], [0.0, 0.0, i_dis], [0.0, i_dis*0.86603, i_dis*0.5]]}
            move = {self.element : [[0, 0, 0], [0, 0, 0], [0, 0, 0]]}
        else:
            numbers = {self.element : 2}
            positions_angstrom_lat0 = {self.element : [[0.0, 0.0, 0.0], [0.0, 0.0, i_dis]]}
            move = {self.element : [[0, 0, 0], [0, 0, 0]]}
        obj = Stru(elements=[self.element], positions_angstrom_lat0=positions_angstrom_lat0, lat0=lat0, cell=cell, pps=pps, numbers=numbers, masses=masses, magmoms=magmoms, move=move)
        obj.write_stru(folder/"STRU")

        return str(folder)

    def _execute(self, dst, command, count=1, **kwargs):
        """Execute calculation

        :params dst: path of working directory
        :params command: string to execute a scf calculation
        :params count: how many times does it take for all dimers to be calculated, each time `(number of dimers)/count` be calculated
        """

        def set_commands(one_list):
            commands = []
            for folder in one_list:
                commands.append(f"(cd {folder};{command};cd ../) &")
            line =  '\n'.join(commands)+"\nwait\n"
            return line
        current_path = Path.cwd()
        lines = []
        all_list = np.array_split(self.folder_list, count)
        for one_list in all_list:
            lines.append(set_commands(one_list))
        os.system('\n'.join(lines))
        os.chdir(current_path)

    def _check(self, dst, **kwargs):
        """Check if job is finished
        
        :params dst: path of working directory
        """
        
        for folder in self.folder_list:
            if not glob.glob(f"{folder}/test.*.dat"):
                raise FileNotFoundError(f"'test.*.dat' file not found in {folder}, may be something wrong.")

    def _parse(self, dst, **kwargs):
        """Parse output of a finished job and write a file named 'folder' for optimize ABFS orbitals
        
        :params dst: path of working directory
        :return: some result
        """

        filename = Path(dst, "folders")
        origin = [str(folder/"test.0.dat") for folder in self.folder_list]
        if self.target == 1:
            linear = [[str(folder/"test.1.dat") for folder in self.folder_list]]
        data = {"origin" : origin, "linear" : linear}
        with open(filename, 'w') as file:  
            json.dump(data, file, indent=4)
        
        res = {"dimer_matrix":1}

        return res

class OptLCAO(ABACUSCalculation):
    """Optimize atomic orbitals(LCAO)"""

    def __init__(self, input_dict, src, element, Nu, rcut, ref_band, target=1, dr=0.01, lr=0.01, **kwargs):
        """Set input parameters of LCAOs optimization
        
        :params input_dict: dict of input parameters
        :params src: path of example which will be tested
        :params element: string of element name
        :params Nu: size of LCAO, e.g. `[4,3,2,1]` means 4s3p2d1f
        :params rcut: cutoff of wavefunctions(a.u.)
        :params ref_band: number of  reference bands
        :params target: fitting target, 0 - wave function, 1 - wave function and its gradient
        """

        super().__init__(input_dict, src, **kwargs)
        if "ecutwfc" not in self.input_dict.keys():
            raise KeyError("`ecutwfc` must be set in `input_dict`.")
        self.element = element
        self.Nu = Nu
        self.rcut = rcut
        self.dr = dr
        self.lr = lr
        self.ref_band = ref_band
        self.target = target
        self.dimer_num = len(distance[self.element])

    def _prepare(self, dst, **kwargs):
        """Prepare input files for optimizing ABFs e.g. input.json
        
        :params dst: path of working directory
        """

        self.folder_opt = Path(dst, "opt_orb_"+"-".join(list_elem2str(self.Nu)))
        self.folder_opt.mkdir(exist_ok=False)

        if not Path(self.src, "folders").exists():
            raise FileNotFoundError("'folders' which is a out file of `SetDimers` calculations not found.")

        with open(self.folder_opt/"input.json", 'w') as file:
            json.dump(self.set_input(), file, indent=4)

    def set_input(self):
        input_params = {
            "file_list": read_json(Path(self.src, "folders")),
            "info":
            {
                "Nt_all":   [self.element],
                "Nu":       {self.element:self.Nu},
                "Nb_true":	[self.ref_band]*self.dimer_num,
                "weight":	[1]*self.dimer_num,
                "Rcut":		{self.element:self.rcut},
                "dr":		{self.element:self.dr},
                "Ecut":		{self.element:self.input_dict["ecutwfx"]},
                "lr":		self.lr
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

    def _execute(self, dst, command, **kwargs):
        """Execute calculation

        :params dst: path of working directory
        """
        
        current_path = Path.cwd()
        os.chdir(self.folder_opt)
        os.environ["MKL_THREADING_LAYER"] = "GNU"
        os.system(command)
        os.chdir(current_path)

    def _check(self, dst, **kwargs):
        """Check if job is finished
        
        :params dst: path of working directory
        """

        filename = f"{self.folder_opt}/orb_{self.element}.dat"
        if not glob.glob(filename):
            raise FileNotFoundError(f"{filename} not found.")

    def _parse(self, dst, **kwargs):
        """Parse output of a finished job
        
        :params dst: path of working directory
        """

        # Orbitals
        filename = f"{self.folder_opt}/orb_{self.element}.dat"
        if self.target == 0:
            datafile = f"{self.folder_opt}/psi_{self.element}.dat"
            name = "PSI"
        elif self.target == 1:
            datafile = f"{self.folder_opt}/psi_{self.element}.dat"
            name = "DPSI"
        os.rename(filename, datafile)

        res = {name:1}
        return res

class LCAO(ABACUSCalculation):
    """Linear Combination of Atomic Orbitals"""

    def __init__(self, input_dict, src, element, Nu, rcut, ref_band, target=1, dr=0.01, lr=0.01, **kwargs) -> None:
        """Set input parameters of dimers calculation

        :params input_dict: dict of input parameters
        :params src: path of example which will be tested
        :params element: string of element name
        :params Nu: size of LCAO, e.g. `[4,3,2,1]` means 4s3p2d1f
        :params rcut: cutoff of wavefunctions(a.u.)
        :params ref_band: number of  reference bands
        :params target: fitting target, 0 - wave function, 1 - wave function and its gradient
        """

        super().__init__(input_dict, src, **kwargs)
        self.element = element
        self.Nu = Nu
        self.rcut = rcut
        self.ref_band = ref_band
        self.target = target
        self.dr = dr
        self.lr = lr
        self.kwargs = kwargs
    
    def _prepare(self, dst, **kwargs):
        """Prepare input files for hybrid  e.g. input.json
        
        :params dst: path of working directory
        """

        self.obj_setdimers = SetDimers(self.input_dict, self.src, self.element, self.Nu, self.rcut, self.ref_band, self.target, **self.kwargs)
        self.obj_optlcao = OptLCAO(self.input_dict, dst, self.element, self.Nu, self.rcut, self.ref_band, self.target, self.dr, self.lr, **self.kwargs)

    def _execute(self, dst, command, **kwargs):
        """Execute calculation

        :params dst: path of working directory
        :params command: pyautotest.schedulers.data.Code` object to execute calculation
        """

        # set dimer
        if isinstance(command, Code):
            dimer_line = Code(code_name=command.code_name, 
                            cmdline_params=command.cmdline_params, 
                            stdout_name=command.stdout_name,
                            stderr_name=command.stderr_name,
                            withmpi=command.withmpi).run_line()
        elif isinstance(command, str):
            dimer_line = command
        self.obj_setdimers.calculate(dst, dimer_line, **kwargs)

        # optimize LCAO
        external_command = kwargs.pop("external_command")
        self.obj_optlcao.calculate(dst, external_command, **kwargs)

    def _check(self, dst, **kwargs):
        """Check if job is finished
        
        :params dst: path of working directory
        """
        pass

    def _parse(self, dst, **kwargs):
        """Parse output of a finished job
        
        :params dst: path of working directory
        :return: some result
        """

        res = {"LCAO":1}
        return res