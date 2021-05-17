'''
Date: 2021-04-02 19:45:22
LastEditors: jiyuyang
LastEditTime: 2021-04-28 18:03:09
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.schedulers.data import Code
from pyautotest.utils.typings import *
from pyautotest.calculations.structure import *
from pyautotest.utils.IO import read_stru, read_kpt

import re
import os
import json
from pathlib import Path
from collections import OrderedDict

def set_cal(name: str, input_dict: dict, stru: Stru, kpt: Kpt, **kwargs):
    module, cal = name.split('.')
    exec(f"""from pyautotest.calculations.plugins.{module} import {cal}""")
    return eval(f"{cal}(input_dict=input_dict, stru=stru, kpt=kpt, **{kwargs})")

def configure(config_file: str_PathLike, version: list) -> Return_2:
    """Configure commands and workflow for Autotest

    :params config_file: configure file for Autotest
    :params version: list of absolute path of executable file
    """
    
    workflow = []
    allcommands = []
    with open(config_file, 'r') as file:
        text = json.load(file)
    src = Path(config_file).parent
    workflowinput = text["workflow"]
    for index, cal in enumerate(workflowinput.keys()):
        if re.match(f"cal_{index}", cal):
            cal_elem = workflowinput[cal]
            name = cal_elem.pop("name", None)
            if not name:
                raise KeyError(f"Key `name` in {cal} isn't set")
            code = cal_elem.pop("code", None)
            cmdline_params = code.pop("cmdline_params", [])
            stdin_name = code.pop("stdin_name", None)
            join_files=code.pop("join_files", False)
            withmpi=code.pop("withmpi", "mpirun")
            commands = []
            for j, ver in enumerate(version):
                commands.append(Code(code_name=ver,
                                    cmdline_params=cmdline_params, 
                                    stdin_name=stdin_name, 
                                    stdout_name=f"{cal}.log", 
                                    stderr_name=f"{cal}.err",
                                    join_files=join_files,
                                    withmpi=withmpi
                                    ))
            allcommands.append(commands)
            input_dict = cal_elem.pop("input_params", None)

            # STRU
            stru_file = src/cal_elem.pop("stru_file", "STRU")
            if stru_file.exists():
                stru = read_stru(input_dict["ntype"], stru_file)
                for elem in stru.elements:
                    if input_dict and "pseudo_dir" not in input_dict.keys():
                        stru.pps[elem] = Path(src, stru.pps[elem])
                    stru.orbitals[elem] = Path(src, stru.orbitals[elem])
            else:
                stru = None

            # KPT
            kpt_file = src/cal_elem.pop("kpt_file", "KPT")
            if kpt_file.exists():
                kpt = read_kpt(kpt_file)
            else:
                kpt = None

            obj = set_cal(name=name, input_dict=input_dict, stru=stru, kpt=kpt, **cal_elem)
            workflow.append(obj)
        else:
            raise KeyError("Wrong keyword `cal` settings, it should be set in numerical order.")
            
    return allcommands, workflow

class Autotest:
    """Auto-test for ABACUS"""

    def __init__(self, workflow: list):
        """Initialize for auto-test

        :params workflow: list, list of class, e.g. [SCF(**kwargs), Band(**kwargs)]
        """

        self.workflow = workflow

    def calculate(self, command: Command, external_command: Command="", save_files: bool=False) -> dict:
        """The whole process of auto-test
        
        :params command: string of command line or `pyautotest.schedulers.data.Code` object.
        :params external_command: other non-ABACUS code needed. Default: ""
        :params save_files: save input files of last calculation or not. Default: False
        """

        for index, cal in enumerate(self.workflow):
            if save_files:
                save_dir = f"cal_{str(index)}"
            else:
                save_dir = ""
            print(f"Begin {cal.__str__()}", flush=True)
            res = cal.calculate(command=command, index=index, external_command=external_command, save_dir=save_dir)
            print(f"End {cal.__str__()}", flush=True)
        return res

    def compare(self, commands: List_Command, external_command: Command="", save_files: bool=False):
        """Comparison test between different commands
        
        :params commands: list of commands string or `pyautotest.schedulers.data.Code` object
        :params external_command: other non-ABACUS code needed Default: ""
        :params save_files: save input files of last calculation or not. Default: False
        """

        res = OrderedDict()
        for index, cal in enumerate(self.workflow):
            print(f"Begin {cal.__str__()}", flush=True)
            for j, command in enumerate(commands[index]):
                subdst = Path(f"command_{j}")
                if not subdst.exists():
                    subdst.mkdir(parents=True)
                os.chdir(subdst)
                if save_files:
                    save_dir = f"cal_{index}"
                else:
                    save_dir = ""
                res[f"command_{j}"] = cal.calculate(command=command, index=index, external_command=external_command, save_dir=save_dir)
                os.chdir("../")
            print(f"End {cal.__str__()}", flush=True)
            self._check(res)
    
    def _check(self, res: dict):
        for index, version in enumerate(res.items()):
            value_version = version[1]
            if index >= 1:
                for key in value_version.keys():
                    if abs(res[f"command_{index}"][key]-res[f"command_{index-1}"][key]) > 1e-8:
                        raise Exception(f"The `{key}` of Version_{index-1} and Version_{index} are different.")
                    else:
                        print(f"The `{key}` of Version_{index-1} and Version_{index} are same.", flush=True)
            else:
                print(f"Only one version calculated.", flush=True)
