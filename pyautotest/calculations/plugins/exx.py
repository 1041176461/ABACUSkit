'''
Date: 2021-03-24 16:05:38
LastEditors: jiyuyang
LastEditTime: 2021-04-29 15:24:09
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.utils.script import set_scheduler
from pyautotest.calculations.structure import *
from pyautotest.calculations.baseclass import ABACUSCalculation, JobCalculation
from pyautotest.calculations.plugins.scf import SCF
from pyautotest.schedulers.data import Code
from pyautotest.utils.tools import *
from pyautotest.utils.typings import *

import re
import os
import glob
import json
import shutil
import typing
import numpy as np
from copy import deepcopy
from pathlib import Path
from collections import defaultdict

class SetDimers(ABACUSCalculation):
    """Set and calculate dimers for off-site Auxiliary Basis Functions(ABFs)"""

    def __init__(self, input_dict: dict, stru: typing.Optional[Stru], kpt: typing.Optional[Kpt], Nu: list, dimer_num: int, **kwargs) -> None:
        """Set input parameters of dimers calculation

        :params input_dict: dict of input parameters
        :params stru: object of `pyautotest.calculations.structure.Stru`
        :params kpt: object of `pyautotest.calculations.structure.Kpt`
        :params Nu: size of ABFs, e.g. `[4,3,2,1]` means 4s3p2d1f
        :params dimer_num: number of dimers
        """

        super().__init__(input_dict, stru, kpt, **kwargs)
        if self.input_dict["basis_type"] not in ["lcao", "lcao_in_pw"]:
            raise KeyError("ABFs based on lcao basis, so `basis_type` in 'INPUT' must be 'lcao' or 'lcao_in_pw'.")
        self.Nu = Nu
        self.dimer_num = dimer_num
        self.folder_weight = {}
        self.orb_obj_list = []
        if  len(self.stru.orbitals) != 0:
            for elem in self.stru.orbitals:
                filename = self.stru.orbitals[elem]
                orb_obj = read_orb(filename)
                self.orb_obj_list.append(orb_obj)
        else: 
            raise FileNotFoundError("Can not find orbital files.")

    def _prepare(self, **kwargs):
        """Prepare input files for dimers calculation"""

        super()._prepare(**kwargs)

        # KPT
        if "gamma_only" not in self.input_dict.keys() or self.input_dict["gamma_only"] == 0:
                kpoints = self.kpt.numbers 
        elif self.input_dict["gamma_only"] == 1:
            kpoints = [1, 1, 1]
        else:
            raise FileNotFoundError("`gamma_only` can only be 1 or 0.")

        # Orbital
        nw_dict = {}
        rcut_dict = {}
        for orb_obj in self.orb_obj_list:
            elem = orb_obj.element
            nw_dict[elem] = orb_obj.total
            rcut_dict[elem] = orb_obj.rcut
        dis = cut_dis(cal_dis(self.stru.supercell_positions(kpoints)), rcut_dict)
        dis_decimal = round_dis(dis, 1e-6)
        dis_opt = self.get_dis_opt(dis_decimal)
        dis_weight = self.cal_dis_weight(dis_opt, dis)
        for T1, T2 in dis_weight:
            for i_dis, i_weight in dis_weight[T1,T2].items():
                self.folder_weight[self.create_input(T1, T2, i_dis, nw_dict)] = i_weight

    def create_input(self, T1: str, T2: str, i_dis: float, nw_list: dict) -> str:
        """Create dimer directory and its input files
        
        :params T1: type 1
        :params T2: type 2
        :params i_dis: distance between `T1` and `T2`
        :params nw_list: dict of total number of orbitals
        """

        folder = Path(folder_name(T1, T2, i_dis))
        folder.mkdir(parents=True, exist_ok=False)
        elements = [T1, ] if T1==T2 else [T1, T2]
        
        # Orbital
        for elem in elements:
            filename = self.stru.orbitals[elem]
            shutil.copy(filename, folder)

        # Pseudopotential
        if "pseudo_dir" not in self.input_dict.keys():
            for elem in elements:
                filename = self.stru.pps[elem]
                shutil.copy(filename, folder)

        input_copy = deepcopy(self.input_dict)
        input_copy["ntype"] = 1 if T1==T2 else 2
        input_copy["calculation"] = "scf"
        input_copy["exx_hybrid_type"] = 'opt_orb'
        input_copy["nbands"] = (nw_list[T1] if abs(i_dis)<1e-10 else nw_list[T1]+nw_list[T2])
        input_copy["nspin"] = 1
        input_copy["gamma_only"] = 1
        input_copy["exx_opt_orb_lmax"] = len(self.Nu)-1

        # INPUT
        input_lines = self.get_input_line(input_copy)
        with open(folder/"INPUT", 'w') as file:
            file.write(input_lines)

        #STRU
        lat0=1
        cell=[[30, 0, 0], [0, 30, 0], [0, 0, 30]]
        if T1 == T2:
            pps = {T1 : Path(self.stru.pps[T1]).name}
            orbitals = {T1 : Path(self.stru.orbitals[T1]).name}
            masses = {T1 : 1}
            magmoms = {T1 : 0}
            if abs(i_dis) < 1e-10:
                positions = {T1 : [[0.0, 0.0, 0.0]]}
                move = {T1 : [[0, 0, 0]]}
            else:
                positions = {T1 : [[0.0, 0.0, 0.0], [i_dis, 0.0, 0.0]]}
                move = {T1 : [[0, 0, 0], [0, 0, 0]]}
        else:
            pps = {T1 : Path(self.stru.pps[T1]).name, T2 : self.stru.pps[T2]}
            orbitals = {T1 : Path(self.stru.orbitals[T1]).name, T2 : Path(self.stru.orbitals[T2]).name}
            masses = {T1 : 1, T2 : 1}
            magmoms = {T1 : 0, T2 : 0}
            positions = {T1 : [[0.0, 0.0, 0.0]], T2 : [[i_dis, 0.0, 0.0]]}
            move = {T1 : [[0, 0, 0]], T2 : [[0, 0, 0]]}
        obj = Stru(positions=positions, lat0=lat0, cell=cell, pps=pps, orbitals=orbitals, masses=masses, magmoms=magmoms, move=move)
        obj.write_stru(folder/"STRU")

        return str(folder)

    def get_dis_opt(self, dis: Dict_Tuple_Dict, opt_mode: str="kmeans") -> Dict_Tuple_Dict:
        """Select some representative point in dictionary of distance
        
        :params dis: dict, distance between two atoms. Its format is `dis[T1,T2] = {..., i_dis:num, ...}`
        :params opt_mode: str, way to select
        """

        dis_opt = dict()
        for T1,T2 in dis:
            dis_tmp = delete_zero(dis[T1,T2])
            if len(dis_tmp)<=self.dimer_num:
                dis_opt[T1,T2] = list(dis_tmp.keys())
            else:
                if opt_mode=="linspace":
                    dis_opt[T1,T2] = np.linspace( min(dis_tmp), max(dis_tmp), self.dimer_num )
                elif opt_mode=="kmeans":
                    from sklearn.cluster import KMeans
                    kmeans = KMeans(n_clusters=self.dimer_num)
                    kmeans.fit_predict(np.array(list(dis_tmp.keys())).reshape(-1,1), 
                                        sample_weight = [num/i_dis**2 for i_dis, num in dis_tmp.items()])
                    dis_opt[T1,T2] = list(kmeans.cluster_centers_.reshape(-1))
            if T1==T2:
                dis_opt[T1,T2].append(0.0)
                
        return dis_opt

    # dis_weight[T1,T2] = {..., i_dis_opt:weight, ...}
    def cal_dis_weight(self, dis_opt: Dict_Tuple_Dict, dis_all: Dict_Tuple_Dict) -> Dict_Tuple_Dict:
        def weight_func(x,D):
            return (2/D**3)*x**3 + (-3/D**2)*x**2 + 1
        dis_weight = defaultdict(dict)
        for T1,T2 in dis_opt:
            dis_opt_TT = sorted(dis_opt[T1,T2])
            for index,i_dis_opt in enumerate(dis_opt_TT):
                i_weight = 0
                i_dis_low = dis_opt_TT[index-1] if index>0 else -np.infty
                i_dis_up = dis_opt_TT[index+1] if index<len(dis_opt_TT)-1 else np.infty
                for i_dis,num in dis_all[T1,T2].items():
                    if i_dis_low<i_dis<i_dis_opt:
                        i_weight += weight_func( i_dis_opt-i_dis, i_dis_opt-i_dis_low ) * num
                    elif i_dis==i_dis_opt:
                        i_weight += num
                    elif i_dis_opt<i_dis<i_dis_up:
                        i_weight += weight_func( i_dis_opt-i_dis, i_dis_opt-i_dis_up ) * num
                dis_weight[T1,T2][i_dis_opt] = i_weight

        return dis_weight

    def _execute(self, command: Command, count: int=1, **kwargs):
        """Execute calculation

        :params command: string to execute a scf calculation
        :params count: how many times does it take for all dimers to be calculated, each time `(number of dimers)/count` be calculated
        """

        filename = "folders"
        #TODO: for new dpsi scripts, it should be a dict with key "origin"
        subfolder_weight = {os.path.basename(key):value for key,value in self.folder_weight.items()}
        with open(filename, 'w') as file:
            json.dump(subfolder_weight, file, indent=4)

        def set_commands(one_list, command):
            if isinstance(command, Code):
                command = command.run_line()
            commands = []
            for folder in one_list:
                commands.append(f"(cd {folder};{command};cd ../) &")
            newline =  '\n'.join(commands)+"\nwait\n"
            return newline
        current_path = Path.cwd()
        lines = []
        all_list = np.array_split(list(self.folder_weight.keys()), count)
        for one_list in all_list:
            lines.append(set_commands(one_list, command))

        os.system('\n'.join(lines))
        os.chdir(current_path)

        # double check
        newcommands = []
        for folder in self.folder_weight:
            if not glob.glob(f"{folder}/matrix_*"):
                import warnings
                warnings.warn("'matrix_*' file not found in {folder}, it will execute again.")
                newcommands.append(f"(cd {folder};{command};cd ../) &")
        os.system('\n'.join(newcommands))
        os.chdir(current_path)

    def batch_run(self, command: Code, scheduler: str, **kwargs):
        """Batch execute dimer calculation
        
        :params command: pyautotest.schedulers.data.Code` object to execute calculation
        :params scheduler: string of scheduler name
        :params kwargs: other parameters of scheduler
        """
        self._prepare()
        filename = "folders"
        #TODO: for new dpsi scripts, it should be a dict with key "origin"
        subfolder_weight = {os.path.basename(key):value for key,value in self.folder_weight.items()}
        with open(filename, 'w') as file:
            json.dump(subfolder_weight, file, indent=4)
        
        current_path = Path.cwd()
        codes_info = [command]
        for folder in self.folder_weight:
            os.chdir(folder)
            submit_command = set_scheduler(scheduler, codes_info, **kwargs)
            os.system(submit_command)
            os.chdir(current_path)

    def double_run(self, command: Command, scheduler: str, **kwargs):
        """Check if job is finished and execute again
        
        :params command: pyautotest.schedulers.data.Code` object to execute calculation
        :params scheduler: string of scheduler name
        :params kwargs: other parameters of scheduler
        """

        current_path = Path.cwd()
        codes_info = [command]
        for folder in read_json("folders"):
            if not glob.glob(f"{folder}/matrix_*"):
                import warnings
                warnings.warn("'matrix_*' file not found in {folder}, it will execute again.")
                os.chdir(folder)
                submit_command = set_scheduler(scheduler, codes_info, **kwargs)
                os.system(submit_command)
                os.chdir(current_path)

    def _check(self, **kwargs):
        """Check if job is finished"""

        for folder in self.folder_weight:
            if not glob.glob(f"{folder}/matrix_*"):
                raise FileNotFoundError(f"'matrix_*' file not found in {folder}, may be something wrong.")

    def _parse(self, **kwargs) -> dict:
        """Parse output of a finished job and write a file named 'folder' for optimize ABFS orbitals
        
        :return: some result
        """

        res = {"dimer_matrix":1}

        return res

class OptABFs(JobCalculation):
    """Optimize off-site Auxiliary Basis Functions(ABFs)"""

    def __init__(self, ecut: float, stru: typing.Optional[Stru], Nu: list, dr: float=0.01, lr: float=0.01, **kwargs):
        """Set input parameters of off-site ABFs optimization
        
        :params ecut: energy cutoff
        :params stru: object of `pyautotest.calculations.structure.Stru`
        :params Nu: size of ABFs, e.g. `[4,3,2,1]` means 4s3p2d1f
        """

        super().__init__(**kwargs)
        self.ecut = ecut
        self.stru = stru
        self.Nu = Nu
        self.dr = dr
        self.lr = lr
        self.orb_obj_list = []
        if  len(self.stru.orbitals) != 0:
            for elem in self.stru.orbitals:
                filename = self.stru.orbitals[elem]
                orb_obj = read_orb(filename)
                self.orb_obj_list.append(orb_obj)
        else: 
            raise FileNotFoundError("Can not find orbital files.")

    def _prepare(self, **kwargs):
        """Prepare input files for optimizing ABFs e.g. input.json"""

        self.folder_opt = Path("opt_orb_"+"-".join(list_elem2str(self.Nu)))
        self.folder_opt_matrix = self.folder_opt/"matrix"
        self.folder_opt_matrix.mkdir(parents=True, exist_ok=False)

        if not Path("folders").exists():
            raise FileNotFoundError("'folders' which is a out file of `SetDimers` calculations not found.")

        self.folders_weight = read_json("folders")
        self.set_matrix()
        with open(self.folder_opt/"input.json", 'w') as file:
            json.dump(self.set_input(), file, indent=4)

    def set_matrix(self):
        for folder_matrix in self.folders_weight:
            T1,T2 = folder_matrix.split("_")[0].split("-")
            dis = float(folder_matrix.split("_")[1])
            if dis:
                if T1==T2:
                    matrix_file = "matrix_0_0_0_1" 
                else:
                    matrix_file = "matrix_0_0_1_0"
            else:
                matrix_file = "matrix_0_0_0_0"
            copy_file = self.folder_opt_matrix/folder_matrix
            shutil.copyfile(os.path.join(folder_matrix, matrix_file), copy_file)
            with open(copy_file, 'r') as file:
                nband = int(re.compile(r"(\d+)\s+nbands").search(file.read()).group(1))
                self.folders_weight[folder_matrix] = (self.folders_weight[folder_matrix], nband)
    
    def set_input(self):
        input_params = {
            "file_list": [str(self.folder_opt_matrix/folder_matrix) for folder_matrix in self.folders_weight], #TODO: for new dpsi scripts, it should be a dict with key "origin"
            "info":
            {
                "Nt_all":       self.stru.elements,
                "Nu":           {elem:self.Nu for elem in self.stru.elements},
                "Nb_true":	    [nbands for weight,nbands in self.folders_weight.values()],
                "weight":	    [weight for weight,nbands in self.folders_weight.values()],
                "Rcut":		    {obj.element:obj.rcut for obj in self.orb_obj_list},
                "dr":		    {elem:self.dr for elem in self.stru.elements},
                "Ecut":		    {elem:self.ecut for elem in self.stru.elements},
                "lr":		    self.lr,
                "cal_T":        "false",
                "cal_smooth":   "false"
            },
            "C_init_info":
            { 
                "init_from_file":   False 
            },
            "V_info":
            {
                "init_from_file":	True,
                "same_band":		False
            }
        }
        return input_params

    def _execute(self, command: Command, **kwargs):
        """Execute calculation"""
        
        if isinstance(command, Code):
            command = command.run_line()

        current_path = Path.cwd()
        os.chdir(self.folder_opt)
        os.environ["MKL_THREADING_LAYER"] = "GNU"
        os.system(command)
        os.chdir(current_path)

    def _check(self, **kwargs):
        """Check if job is finished"""

        for elem in self.stru.elements:
            filename = f"{self.folder_opt}/orb_{elem}.dat"
            if not glob.glob(filename):
                raise FileNotFoundError(f"{filename} not found.")

    def _parse(self, **kwargs):
        """Parse output of a finished job"""

        subdst = Path("exx_"+"-".join(list_elem2str(self.Nu)))
        subdst.mkdir(parents=True, exist_ok=False)

        #ABFs
        for elem in self.stru.elements:
            filename = f"{self.folder_opt}/orb_{elem}.dat"
            datafile = f"{subdst}/abfs_{elem}.dat"
            self.stru.abfs[elem] = datafile
            shutil.copyfile(filename, datafile)

        # STRU
        self.stru.write_stru(f"{subdst}/STRU")
        
        # Orbitals
        for obj in self.orb_obj_list:
            shutil.copy(obj.datafile, subdst)

        # Pseudopotential
        if self.stru.pps:
            for i in self.stru.pps.values():
                shutil.copy(i, subdst)

        res = {"ABFs":1}
        return res

#TODO: Now EXX only support scf calculations, so we inherit class SCF, maybe inherit RELAX or CELL_RELAX later on 
class EXX(SCF):
    """Hybrid function calculation with off-site ABFs basis"""

    def __init__(self, input_dict: dict, stru: typing.Optional[Stru], kpt: typing.Optional[Kpt], Nu: list, dimer_num: int, dr: float=0.01, lr: float=0.01, **kwargs) -> None:
        """Set input parameters of hybrid function calculation

        :params input_dict: dict of input parameters
        :params stru: object of `pyautotest.calculations.structure.Stru`
        :params kpt: object of `pyautotest.calculations.structure.Kpt`
        :params Nu: size of ABFs, e.g. `[4,3,2,1]` means 4s3p2d1f
        :params dimer_num: number of dimers
        """

        super().__init__(input_dict, stru, kpt, **kwargs)
        if self.input_dict["basis_type"] not in ["lcao", "lcao_in_pw"]:
            raise KeyError("ABFs based on lcao basis, so `basis_type` in 'INPUT' must be 'lcao' or 'lcao_in_pw'.")
        self.delete_key()
        self.Nu = Nu
        self.dimer_num = dimer_num
        self.dr = dr
        self.lr = lr
        self.kwargs = kwargs

    def delete_key(self):
        key_list = ["ocp", "ocp_set", "nelec"]
        for i in key_list:
            self.input_dict.pop(i, None)

    def _prepare(self, **kwargs):
        """Prepare input files for hybrid  e.g. input.json"""

        self.obj_setdimers = SetDimers(self.input_dict, self.stru, self.kpt, self.Nu, self.dimer_num, **self.kwargs)
        self.obj_optabfs = OptABFs(self.input_dict["exx_opt_orb_ecut"], self.stru, self.Nu, self.dr, self.lr, **self.kwargs)

    def _execute(self, command: Code, external_command: Command, **kwargs):
        """Execute calculation

        :params command: pyautotest.schedulers.data.Code` object to execute calculation
        """
        
        if not isinstance(command, Code):
            raise TypeError("Type of `command` here must be `pyautotest.schedulers.data.Code` object")

        current_path = Path.cwd()
        # set dimer
        dimer_line = Code(code_name=command.code_name, 
                        cmdline_params=['-np 1'], 
                        stdout_name=command.stdout_name,
                        stderr_name=command.stderr_name,
                        withmpi=command.withmpi).run_line()
        self.obj_setdimers.calculate(dimer_line, **kwargs)

        # optimize ABFs
        self.obj_optabfs.calculate(external_command, **kwargs)
        
        # Exx calculation
        os.chdir("exx_"+"-".join(list_elem2str(self.Nu)))
        # INPUT
        input_lines = self._get_input_line()
        with open("INPUT", 'w') as file:
            file.write(input_lines)
        os.system(command.run_line())
        os.chdir(current_path)
