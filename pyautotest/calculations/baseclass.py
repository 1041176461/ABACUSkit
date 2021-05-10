'''
Date: 2021-03-05 17:03:37
LastEditors: jiyuyang
LastEditTime: 2021-04-28 15:50:39
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.schedulers.data import Code
from pyautotest.calculations.structure import Stru, Kpt
from pyautotest.utils.typings import *

import os
import re
import abc
import shutil
import typing
from pathlib import Path

class JobCalculation(abc.ABC):
    """Single job calculation"""

    def __init__(self, **kwargs) -> None:
        """Set input parameters of single job calcultion"""

        self.kwargs = kwargs

    def __str__(self) -> str:
        return self.__class__.__name__

    @abc.abstractmethod
    def _prepare(self, **kwargs):
        """Prepare input files for a job"""

    def _execute(self, command: Command, **kwargs):
        """Execute calculation

        :params command: string or `pyautotest.schedulers.data.Code` object to execute calculation
        """

        if isinstance(command, Code):
            command = command.run_line()

        os.system(command)

    @abc.abstractmethod
    def _check(self, **kwargs):
        """Check if job is finished"""

    @abc.abstractmethod
    def _parse(self, **kwargs) -> dict:
        """Parse output of a finished job
        
        :return: some result
        """

        res = {}
        return res

    def calculate(self, command: Command, save_dir: str_PathLike="", **kwargs):
        """The whole process of job calculation

        :params command: string or `pyautotest.schedulers.data.Code` object to execute calculation
        """

        self._prepare(**kwargs)
        self._execute(command, **kwargs)
        self._check(**kwargs)
        res = self._parse(**kwargs)
        if save_dir:
            self.save(save_dir)
        return res

    def save(self, save_dir: str_PathLike):
        """Save all input and out files of last calculation
        
        :params save_dir: directory where to save all input and output files
        """

        all_files =  Path(".").iterdir()
        Path(save_dir).mkdir(parents=True, exist_ok=False)
        for file in all_files:
            if Path(file).is_file():
                shutil.copy(file, save_dir)

class ABACUSCalculation(JobCalculation):
    """ABACUS Calculation"""

    def __init__(self, input_dict: dict, stru: typing.Optional[Stru], kpt: typing.Optional[Kpt], **kwargs) -> None:
        """Set input parameters of ABACUS calcultion
        
        :params input_dict: dict of input parameters
        :params stru: object of `pyautotest.calculations.structure.Stru`
        :params kpt: object of `pyautotest.calculations.structure.Kpt`
        """

        super().__init__(**kwargs)
        self.input_dict = input_dict
        self.input_dict["suffix"] = "test"
        if "ntype" not in self.input_dict.keys():
            raise KeyError("Need to set input_dict['ntype']")
        self.stru = stru
        self.kpt = kpt

    def _prepare(self, **kwargs):
        """Prepare input files for ABACUS calculation"""

        # INPUT
        input_lines = self._get_input_line()
        with open("INPUT", 'w') as file:
            file.write(input_lines)

        # STRU
        self.stru.write_stru("STRU")

        # KPT
        if "gamma_only" not in self.input_dict.keys() or self.input_dict["gamma_only"] == 0:
            self.kpt.write_kpt("KPT")
        elif self.input_dict["gamma_only"] != 1:
            raise FileNotFoundError("`gamma_only` can only be 1 or 0.")

    def _check(self, index:int=0, **kwargs) -> typing.Union[int, str]:
        """Check if job is finished
        
        :params index: calculation index in workflow
        """

        time = 0
        with open(f"cal_{index}.log", 'r') as file:
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