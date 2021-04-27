'''
Date: 2021-03-24 16:05:38
LastEditors: jiyuyang
LastEditTime: 2021-04-23 18:06:40
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.calculations.structure import Stru, read_stru, read_orb, read_kpt, cal_dis, cut_dis, round_dis, delete_zero
from pyautotest.calculations.baseclass import ABACUSCalculation
from pyautotest.calculations.plugins.scf import SCF
from pyautotest.utils.tools import read_json, folder_name, list_elem2str
from pyautotest.utils.script import configure_code

import re
import os
import glob
import json
import shutil
import numpy as np
from copy import deepcopy
from pathlib import Path
from collections import defaultdict

class SetDimers(ABACUSCalculation):
    """Set and calculate dimers for off-site Auxiliary Basis Functions(ABFs)"""

    def __init__(self, input_dict, src, Nu, dimer_num, **kwargs) -> None:
        """Set input parameters of dimers calculation

        :params input_dict: dict of input parameters
        :params src: path of example which will be tested
        :params Nu: size of ABFs, e.g. `[4,3,2,1]` means 4s3p2d1f
        :params dimer_num: number of dimers
        """

        super().__init__(input_dict, src, **kwargs)
        if self.input_dict["basis_type"] not in ["lcao", "lcao_in_pw"]:
            raise KeyError("ABFs based on lcao basis, so `basis_type` in 'INPUT' must be 'lcao' or 'lcao_in_pw'.")
        self.Nu = Nu
        self.dimer_num = dimer_num
        self.folder_weight = {}
        self.stru_obj = read_stru(self.input_dict["ntype"], Path(self.src, self.stru_file))
        self.orb_obj_list = []
        if  len(self.orbitals) != 0:
            for elem in self.orbitals:
                filename = Path(self.src, self.orbitals[elem])
                orb_obj = read_orb(filename)
                self.orb_obj_list.append(orb_obj)
        else: 
            raise FileNotFoundError("Can not find orbital files.")

    def _prepare(self, dst, **kwargs):
        """Prepare input files for dimers calculation e.g. INPUT, STRU, KPT, orbital and pseudopotential files
        
        :params dst: path of working directory
        """

        super()._prepare(dst, **kwargs)

        # KPT
        if "gamma_only" not in self.input_dict.keys() or self.input_dict["gamma_only"] == 0:
            if "kpoint_file" in self.input_dict.keys():
                kpt, offset = read_kpt(Path(self.src, self.input_dict["kpoint_file"]))
            else:
                kpt, offset = read_kpt(Path(self.src, self.kpt_file))
        elif self.input_dict["gamma_only"] == 1:
            kpt, offset = [1.0, 1.0, 1.0], [0, 0, 0]
        else:
            raise FileNotFoundError("Can not find k-points file.")

        # Orbital
        nw_dict = {}
        rcut_dict = {}
        for orb_obj in self.orb_obj_list:
            elem = orb_obj.element
            nw_dict[elem] = orb_obj.total
            rcut_dict[elem] = orb_obj.rcut
        dis = cut_dis(cal_dis(self.stru_obj.supercell_positions(kpt)), rcut_dict)
        dis_decimal = round_dis(dis, 1e-6)
        dis_opt = self.get_dis_opt(dis_decimal)
        dis_weight = self.cal_dis_weight(dis_opt, dis)
        for T1, T2 in dis_weight:
            for i_dis, i_weight in dis_weight[T1,T2].items():
                self.folder_weight[self.create_input(T1, T2, i_dis, nw_dict, dst)] = i_weight

    def create_input(self, T1, T2, i_dis, nw_list, dst):
        """Create dimer directory and its input files
        
        :params T1: type 1
        :params T2: type 2
        :params i_dis: distance between `T1` and `T2`
        :params nw_list: dict of total number of orbitals
        :params dst: directory where to create subdirectory named `T1`-`T2`_`i_dis`
        """

        folder = Path(dst, folder_name(T1, T2, i_dis))
        folder.mkdir(parents=True, exist_ok=False)
        elements = [T1, ] if T1==T2 else [T1, T2]
        
        # Orbital
        for elem in elements:
            filename = Path(self.src, self.orbitals[elem])
            shutil.copyfile(filename, Path(folder, self.orbitals[elem]))

        # Pseudopotential
        if len(self.pps) != 0:
            for elem in elements:
                filename = Path(self.src, self.pps[elem])
                shutil.copyfile(filename, Path(folder, self.pps[elem]))
        elif not self.input_dict["pseudo_dir"]:
            raise FileNotFoundError(f"Can not find pseudopotential files in {self.src}.")

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
        numbers = {}
        positions = {}
        move = {}
        if T1 == T2:
            pps = {T1 : self.pps[T1]}
            orbitals = {T1 : self.orbitals[T1]}
            masses = {T1 : 1}
            magmoms = {T1 : 0}
            if abs(i_dis) < 1e-10:
                numbers = {T1 : 1}
                positions = {T1 : [[0.0, 0.0, 0.0]]}
                move = {T1 : [[0, 0, 0]]}
            else:
                numbers = {T1 : 2}
                positions = {T1 : [[0.0, 0.0, 0.0], [i_dis, 0.0, 0.0]]}
                move = {T1 : [[0, 0, 0], [0, 0, 0]]}
        else:
            pps = {T1 : self.pps[T1], T2 : self.pps[T2]}
            orbitals = {T1 : self.orbitals[T1], T2 : self.orbitals[T2]}
            masses = {T1 : 1, T2 : 1}
            magmoms = {T1 : 0, T2 : 0}
        obj = Stru(elements=elements, positions=positions, lat0=lat0, cell=cell, pps=pps, orbitals=orbitals, numbers=numbers, masses=masses, magmoms=magmoms, move=move)
        obj.write_stru(folder/"STRU")

        return str(folder)

    def get_dis_opt(self, dis, opt_mode="kmeans"):
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
    def cal_dis_weight(self, dis_opt, dis_all):
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

    def _execute(self, dst, command, count=1, **kwargs):
        """Execute calculation

        :params dst: path of working directory
        :params command: string to execute a scf calculation
        :params count: how many times does it take for all dimers to be calculated
        """

        def set_commands(one_list):
            commands = []
            for folder in one_list:
                commands.append(f"(cd {folder};{command};cd ../) &")
            line =  '\n'.join(commands)+"\nwait\n"
            return line
        current_path = Path.cwd()
        lines = []
        all_list = np.array_split(list(self.folder_weight.keys()), count)
        for one_list in all_list:
            lines.append(set_commands(one_list))
        os.system('\n'.join(lines))
        os.chdir(current_path)

    def _check(self, dst, **kwargs):
        """Check if job is finished
        
        :params dst: path of working directory
        """
        
        for folder in self.folder_weight:
            if not glob.glob(f"{folder}/matrix_*"):
                raise FileNotFoundError(f"'matrix_*' file not found in {folder}, may be something wrong.")

    def _parse(self, dst, **kwargs):
        """Parse output of a finished job and write a file named 'folder' for optimize ABFS orbitals
        
        :params dst: path of working directory
        :return: some result
        """

        filename = Path(dst, "folders")
        subfolder_weight = {os.path.basename(key):value for key,value in self.folder_weight.items()}
        with open(filename, 'w') as file:
            json.dump(subfolder_weight, file, indent=4)

        res = {"dimer_matrix":1}
        return res

class OptABFs(ABACUSCalculation):
    """Optimize off-site Auxiliary Basis Functions(ABFs)"""

    def __init__(self, input_dict, src, Nu, mainpy="", dr=0.01, lr=0.01, **kwargs):
        """Set input parameters of  off-site ABFs optimization
        
        :params input_dict: dict of input parameters
        :params src: path of example which will be tested
        :params Nu: size of ABFs, e.g. `[4,3,2,1]` means 4s3p2d1f
        :params mainpy: mainpy to optimize ABFs
        """

        super().__init__(input_dict, src, **kwargs)
        self.Nu = Nu
        self.mainpy = mainpy
        self.stru_obj = read_stru(self.input_dict["ntype"], Path(self.src, self.stru_file))
        
        self.orb_obj_list = []
        if  len(self.orbitals) != 0:
            for elem in self.orbitals:
                filename = Path(self.src, self.orbitals[elem])
                orb_obj = read_orb(filename)
                self.orb_obj_list.append(orb_obj) 
        else: 
            raise FileNotFoundError("Can not find orbital files.")
        self.dr = dr
        self.lr = lr

    def _prepare(self, dst, **kwargs):
        """Prepare input files for optimizing ABFs e.g. input.json
        
        :params dst: path of working directory
        """

        self.folder_opt = Path(dst, "opt_orb_"+"-".join(list_elem2str(self.Nu)))
        self.folder_opt_matrix = self.folder_opt/"matrix"
        self.folder_opt_matrix.mkdir(parents=True, exist_ok=False)

        if Path(self.src, "folders").exists():
            Dir = self.src
        else:
            raise FileNotFoundError("'folders' which is a out file of `SetDimers` calculations not found.")

        self.folders_weight = read_json(Path(Dir, "folders"))
        self.set_matrix(Dir)
        with open(self.folder_opt/"input.json", 'w') as file:
            json.dump(self.set_input(), file, indent=4)

    def set_matrix(self, Dir):
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
            shutil.copyfile(Path(Dir, folder_matrix)/matrix_file, copy_file)
            with open(copy_file, 'r') as file:
                nband = int(re.compile(r"(\d+)\s+nbands").search(file.read()).group(1))
                self.folders_weight[folder_matrix] = (self.folders_weight[folder_matrix], nband)
    
    def set_input(self):
        input_params = {
            "file_list": [str(self.folder_opt_matrix/folder_matrix) for folder_matrix in self.folders_weight],
            "info":
            {
                "Nt_all":   self.stru_obj.elements,
                "Nu":       {elem:self.Nu for elem in self.stru_obj.elements},
                "Nb_true":	[nbands for weight,nbands in self.folders_weight.values()],
                "weight":	[weight for weight,nbands in self.folders_weight.values()],
                "Rcut":		{obj.element:obj.rcut for obj in self.orb_obj_list},
                "dr":		{elem:self.dr for elem in self.stru_obj.elements},
                "Ecut":		{elem:self.input_dict["exx_opt_orb_ecut"] for elem in self.stru_obj.elements},
                "lr":		self.lr
            },
            "C_init_info":{ "init_from_file": False },
            "V_info":
            {
                "init_from_file":	True,
                "same_band":		False
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
        os.system(self.mainpy)
        os.chdir(current_path)

    def _check(self, dst, **kwargs):
        """Check if job is finished
        
        :params dst: path of working directory
        """

        for elem in self.stru_obj.elements:
            filename = f"{self.folder_opt}/orb_{elem}.dat"
            if not glob.glob(filename):
                raise FileNotFoundError(f"{filename} not found.")

    def _parse(self, dst, **kwargs):
        """Parse output of a finished job
        
        :params dst: path of working directory
        """

        subdst = Path(dst, "exx_"+"-".join(list_elem2str(self.Nu)))
        subdst.mkdir(parents=True, exist_ok=False)

        #ABFs
        for elem in self.stru_obj.elements:
            filename = f"{self.folder_opt}/orb_{elem}.dat"
            datafile = f"{subdst}/abfs_{elem}.dat"
            shutil.copyfile(filename, datafile)

        #STRU
        self.stru_obj.write_stru(f"{subdst}/STRU")
        for obj in self.orb_obj_list:
            shutil.copy2(obj.datafile, subdst)

        # KPT
        filename = Path(self.src, self.kpt_file)
        if filename.exists():
            shutil.copyfile(filename, Path(subdst, "KPT"))

        # Orbitals
        if self.orbitals:
            for i in self.orbitals.values():
                shutil.copyfile(Path(self.src, i), Path(subdst, i))

        # Pseudopotential
        if self.pps:
            for i in self.pps.values():
                shutil.copyfile(Path(self.src, i), Path(subdst, i))

        res = {"ABFs":1}
        return res

#TODO: Now EXX only support scf calculations, so we inherit class SCF 
class EXX(SCF):
    """Hybrid function calculation with off-site ABFs basis"""

    def __init__(self, input_dict, src, Nu, dimer_num, mainpy="", dr=0.01, lr=0.01, **kwargs) -> None:
        """Set input parameters of hybrid function calculation

        :params input_dict: dict of input parameters
        :params src: path of example which will be tested
        :params Nu: size of ABFs, e.g. `[4,3,2,1]` means 4s3p2d1f
        :params dimer_num: number of dimers
        """

        super().__init__(input_dict, src, **kwargs)
        if self.input_dict["basis_type"] not in ["lcao", "lcao_in_pw"]:
            raise KeyError("ABFs based on lcao basis, so `basis_type` in 'INPUT' must be 'lcao' or 'lcao_in_pw'.")
        self.Nu = Nu
        self.dimer_num = dimer_num
        self.mainpy = mainpy
        self.dr = dr
        self.lr = lr
        self.kwargs = kwargs
        self.obj_setdimers = SetDimers(input_dict, src, Nu, dimer_num, **kwargs)
        self.obj_optabfs = OptABFs(input_dict, src, Nu, mainpy, dr, lr, **kwargs)

    def _prepare(self, dst, **kwargs):
        """Prepare input files for hybrid  e.g. input.json
        
        :params dst: path of working directory
        """

        super()._prepare(dst, **kwargs)
        self.obj_setdimers = SetDimers(self.input_dict, dst, self.Nu, self.dimer_num, **self.kwargs)
        self.obj_optabfs = OptABFs(self.input_dict, self.src, self.Nu, self.mainpy, self.dr, self.lr, **self.kwargs)

    def _execute(self, dst, command, **kwargs):
        """Execute calculation

        :params dst: path of working directory
        :params command: string to execute a scf calculation
        """

        current_path = Path.cwd()
        self.obj_setdimers.calculate(dst, command, **kwargs) #TODO: how to force dimers calculation with one process
        self.obj_optabfs.calculate(dst, command, **kwargs)
        os.chdir(Path(dst, "exx_"+"-".join(list_elem2str(self.Nu))))
        os.system(command)
        os.chdir(current_path)
