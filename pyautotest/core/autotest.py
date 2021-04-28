'''
Date: 2021-04-02 19:45:22
LastEditors: jiyuyang
LastEditTime: 2021-04-28 18:03:09
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.schedulers.data import Code

import re
import json
from pathlib import Path
from collections import OrderedDict

def set_cal(name="", **kwargs):
    module, cal = name.split('.')
    exec(f"""from pyautotest.calculations.plugins.{module} import {cal}""")
    return eval(f"{cal}(**{kwargs})")

def configure(config_file, src="", dst="", version=[]):
    """Configure parameters for Autotest
    
    :params src: str, path of source directory
    :params dst: str, path of working directory
    :params config_file: configure file for Autotest
    :params version: list of absolute path of executable file
    """
    
    workflow = []
    allcommands = []
    with open(config_file, 'r') as file:
        text = json.load(file)
    workflowinput = text["workflow"]
    independent = text.pop("independent", False)
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
            for ver in version:
                commands.append(Code(code_name=ver,
                                    cmdline_params=cmdline_params, 
                                    stdin_name=stdin_name, 
                                    stdout_name=f"{cal}.log", 
                                    stderr_name=f"{cal}.err",
                                    join_files=join_files,
                                    withmpi=withmpi
                                    ))
            allcommands.append(commands)
            input_dict = cal_elem.pop("input_params")
            if independent:
                obj = set_cal(name, input_dict=input_dict, src=src, **cal_elem)
            else:
                if index == 0:
                    obj = set_cal(name, input_dict=input_dict, src=src, **cal_elem)
                else:
                    obj = set_cal(name, input_dict=input_dict, src=dst, **cal_elem)
            workflow.append(obj)
        else:
            raise KeyError("Wrong keyword `cal` settings, it should be set in numerical order.")
            
    return allcommands, workflow

class Autotest:
    """Auto-test for ABACUS"""

    def __init__(self, dst="", workflow=[]):
        """Initialize for auto-test
        
        :params dst: str, path of working directory
        :params workflow: list, list of class, e.g. [SCF(**kwargs), Band(**kwargs)]
        """
        
        self.dst = dst
        self.workflow = workflow

    def calculate(self, command, save_files=False, external_command=''):
        """The whole process of auto-test
        
        :params command: string of command line or `pyautotest.schedulers.data.Code` object. Default: ""
        :params save_files: save input files of last calculation or not. Default: False
        """

        subdst = Path(self.dst)
        subdst.mkdir(parents=True, exist_ok=False)
        for index, cal in enumerate(self.workflow):
            res = cal.calculate(dst=subdst, command=command, index=index, external_command=external_command)
            if save_files:
                cal.save(dst=subdst, save_dir=Path(subdst, f"cal_{str(index)}"))
                
        return res

    def compare(self, commands=[], save_files=False, external_command=''):
        """Comparison test between different commands
        
        :params commands: list of commands string or `pyautotest.schedulers.data.Code` object
        :params save_files: save input files of last calculation or not. Default: False
        :params external_command: other non-ABACUS code needed
        """

        res = OrderedDict()
        for index, cal in enumerate(self.workflow):
            for j, command in enumerate(commands[index]):
                subdst = Path(self.dst, f"command_{j}")
                subdst.mkdir(parents=True, exist_ok=True)
                res[f"command_{j}"] = cal.calculate(dst=subdst, command=command, index=index, external_command=external_command)
                if save_files:
                    cal.save(dst=subdst, save_dir=Path(subdst, f"cal_{str(index)}"))
            self._check(res)
    
    def _check(self, res):
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
