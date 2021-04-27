'''
Date: 2021-03-29 21:33:09
LastEditors: jiyuyang
LastEditTime: 2021-04-23 16:38:39
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

import os
from pyautotest.utils.tools import read_json

class Show:
    """Show auto-test information"""

    @classmethod
    def show_libinfo(cls, src):
        """Show example library information
    
        :params src: path of library
        """

        print("--------------------------Library Information--------------------------")
        print(f"Path: {src}")
        print("Directory Structure:")
        for index, subsrc in enumerate(os.listdir(src)):
            if os.path.isdir(os.path.join(src, subsrc)) and subsrc != "OUT.test":
                line = f" ({index+1}) " + subsrc
                print(f"{line}")
                subline = '\t' + ', '.join(os.listdir(os.path.join(src, subsrc)))
                print(subline)
            elif os.path.isfile(os.path.join(src, subsrc)):
                subline = '\t' + ', '.join(os.listdir(src))
                print(subline)
                break
            else:
                raise FileNotFoundError(f"No information to show!")

    @classmethod
    def show_cmdline(cls, args):
        if args.lib:
            text = read_json(args.lib)
            cls.show_libinfo(text["src"])
