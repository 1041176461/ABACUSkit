'''
Date: 2021-03-05 17:03:37
LastEditors: jiyuyang
LastEditTime: 2021-04-23 15:50:10
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

import os
import re
import abc
import shutil
from pathlib import Path

class JobCalculation(abc.ABC):
    """Single job calculation"""

    def __init__(self, input_dict, src, **kwargs) -> None:
        """Set input parameters of single job calcultion"""

        self.input_dict = input_dict
        self.src = Path(src)

    @abc.abstractmethod
    def _prepare(self, dst, **kwargs):
        """Prepare input files for a job e.g. INPUT, STRU, KPT, orbital and pseudopotential files
        
        :params dst: path of working directory
        """

    def _execute(self, dst, command, **kwargs):
        """Execute calculation

        :params dst: path of working directory
        :params command: string to execute a scf calculation
        """

        current_path = Path.cwd()
        os.chdir(dst)
        os.system(command)
        os.chdir(current_path)

    @abc.abstractmethod
    def _check(self, dst, **kwargs):
        """Check if job is finished
        
        :params dst: path of working directory
        """

    @abc.abstractmethod
    def _parse(self, dst, **kwargs):
        """Parse output of a finished job
        
        :params dst: path of working directory
        :return: some result
        """

        res = {}
        return res

    def get_input_dict(self):
        """Return a dictionary of input parameters"""
        return self.input_dict

    def calculate(self, dst, command, **kwargs):
        """The whole process of job calculation"""

        self._prepare(dst, **kwargs)
        self._execute(dst, command, **kwargs)
        self._check(dst, **kwargs)
        res = self._parse(dst, **kwargs)
        return res

class ABACUSCalculation(JobCalculation):
    """ABACUS Calculation"""

    def __init__(self, input_dict, src, **kwargs) -> None:
        """Set input parameters of ABACUS calcultion
        
        :params input_dict: dict of input parameters
        :params src: path of example which will be tested
        """

        super().__init__(input_dict, src, **kwargs)
        self.input_dict["suffix"] = "test"
        if "ntype" not in self.input_dict.keys():
            raise KeyError("Need to set input_dict['ntype']")

        if "orbitals" in kwargs.keys():
            self.orbitals = kwargs["orbitals"]
        else:
            self.orbitals = {}

        if "pps" in kwargs.keys():
            self.pps = kwargs["pps"]
        else:
            self.pps = {}
            
        if "stru_file" in kwargs.keys():
            self.stru_file = kwargs["stru_file"]
        else:
            self.stru_file = "STRU"
            
        if "kpt_file" in kwargs.keys():
            self.kpt_file = kwargs["kpt_file"]
        else:
            self.kpt_file = "KPT"

    def _prepare(self, dst, **kwargs):
        """Prepare input files for ABACUS calculation e.g. INPUT, STRU, KPT, orbital and pseudopotential files
        
        :params dst: path of working directory
        """

        # INPUT
        input_lines = self._get_input_line()
        with open(Path(dst, "INPUT"), 'w') as file:
            file.write(input_lines)

        # STRU
        shutil.copyfile(Path(self.src, self.stru_file), Path(dst, "STRU"))

        # KPT
        if "gamma_only" not in self.input_dict.keys() or self.input_dict["gamma_only"] == 0:
            if "kpoint_file" in self.input_dict.keys():
                shutil.copyfile(Path(self.src, self.input_dict["kpoint_file"]), Path(dst, self.input_dict["kpoint_file"]))
            else:
                shutil.copyfile(Path(self.src, self.kpt_file), Path(dst, "KPT"))
        elif self.input_dict["gamma_only"] != 1:
            raise FileNotFoundError("Can not find k-points file.")

        # Orbital
        if  self.orbitals:
            for i in self.orbitals.values():
                shutil.copyfile(Path(self.src, i), Path(dst, i))
        elif self.input_dict["basis_type"] in ["lcao", "lcao_in_pw"]:
            raise FileNotFoundError(f"Can not find orbital files in {self.src}.")
        
        # Pseudopotential
        if self.pps:
            for i in self.pps.values():
                shutil.copyfile(Path(self.src, i), Path(dst, i))
        elif not self.input_dict["pseudo_dir"]:
            raise FileNotFoundError(f"Can not find pseudopotential files in {self.src}.")

    def _check(self, dst, index=0, **kwargs):
        """Check if job is finished
        
        :params dst: path of working directory
        :params index: calculation index in workflow
        """

        time = 0
        with open(Path(dst, f"cal_{index}.log"), 'r') as file:
            for line in file:
                if re.search("TOTAL  Time  : ", line):
                    time = re.search("(TOTAL  Time  : )([a-z0-9\s]+)", line).group(2)
        if time == 0:
            raise Exception(f"Calculation {index} may not finish.")
        
        return time

    def _get_input_line(self):
        """Return input lines in INPUT file"""

        lines = []
        lines.append("INPUT_PARAMETERS")
        for key, value in self.input_dict.items():
            if value:
                lines.append(f"{key.ljust(30)}{value}")
        return '\n'.join(lines)
    
    @staticmethod
    def get_input_line(input_dict):
        """Return input lines in INPUT file
        
        :params input_dict: dict of input parameters
        """

        lines = []
        lines.append("INPUT_PARAMETERS")
        for key, value in input_dict.items():
            if value:
                lines.append(f"{key.ljust(30)}{value}")
        return '\n'.join(lines)

    def save(self, dst, save_dir):
        """Save all input and out files of last calculation
        
        :params dst: path of working directory
        :params save_dir: directory where to save all input and output files
        """

        all_files =  Path(dst).iterdir()
        Path(save_dir).mkdir(parents=True, exist_ok=False)
        for file in all_files:
            if Path(file).is_file():
                shutil.copy2(file, save_dir)