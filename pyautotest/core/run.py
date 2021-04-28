'''
Date: 2021-03-15 15:35:22
LastEditors: jiyuyang
LastEditTime: 2021-04-28 18:57:32
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.core.autotest import Autotest, configure
from pyautotest.utils.script import submit_script
from pyautotest.utils.tools import read_json, write_json

import os

class Run:
    """Ways to run Auto-test code"""

    @classmethod 
    def single_run(cls, src="", dst="", version=[], save_files=True, external_command=''):
        """This function tests an example individually, this example directory should have a configuration file named `config.json`.
        
        :params src: path of library
        :params dst: path of working directory
        :params version: executable file list
        :params save_files: if save input files of last calculation, or not
        :params external_command: other non-ABACUS code needed
        """

        config_file = os.path.join(src, "config.json")
        if os.path.exists(config_file):
            print(f"Test For Example {src} Begin:", flush=True)
            commands, workflow = configure(config_file, src, dst, version)
            job = Autotest(dst, workflow)
            job.compare(commands, save_files, external_command)
            print(f"Test For Example {src} Finished\n", flush=True)
        else:
            raise FileNotFoundError(f"Can't find `config.json` in {src}")

    @classmethod
    def batch_run(cls, src="", dst="", version=[], save_files=True, external_command=''):
        """If there is a library of examples for test, this function can run all the test in a serial way. Each example directory
        should have a configuration file named `config.json`.
        
        :params src: path of library
        :params dst: path of working directory
        :params version: executable file list
        :params save_files: if save input files of last calculation, or not
        :params external_command: other non-ABACUS code needed
        """

        for subsrc in os.listdir(src):
            abs_subsrc= os.path.join(src, subsrc)
            subdst = os.path.join(dst, os.path.basename(subsrc))
            cls.single_run(abs_subsrc, subdst, version, save_files, external_command)

    @classmethod
    def batch_with_script(cls, filename, parallel=False):
        """Batch run with script
        
        :params filename: absolute path of `input.json`
        :params parallel: whether to test in parallel. Default: False
        """

        if parallel:
                text = read_json(filename)
                for subsrc in os.listdir(text["src"]):
                    subdst = os.path.join(text["dst"], os.path.basename(subsrc))
                    os.makedirs(subdst)
                    new_filename = os.path.join(subdst, os.path.basename(filename))
                    write_json(filename, new_filename, src=os.path.join(text["src"], os.path.basename(subsrc)), dst=subdst)
                    os.chdir(subdst)
                    line = f"autotest run --single={new_filename} --local True"
                    submit_script(new_filename, line)

        else:
            line = f"autotest run --batch={filename} --local True"
            submit_script(filename, line)

    @classmethod
    def single_with_script(cls, filename):
        """Batch run with script
        
        :params filename: absolute path of `input.json`
        """

        line = f"autotest run --single={filename} --local True"
        submit_script(filename, line)

    @classmethod
    def run_cmdline(cls, args):
        if args.batch and args.local:
            text = read_json(args.batch)
            save_files = text.pop("save_files", True)
            external_command = text.pop("external_command", '')
            cls.batch_run(text["src"], text["dst"], text["version"], save_files, external_command)

        elif args.batch and not args.local:
            cls.batch_with_script(args.batch, args.parallel)
    
        elif args.single and args.local:
            text = read_json(args.single)
            save_files = text.pop("save_files", True)
            external_command = text.pop("external_command", '')
            cls.single_run(text["src"], text["dst"], text["version"], save_files, external_command)
    
        elif args.single and not args.local:
            cls.single_with_script(args.single)