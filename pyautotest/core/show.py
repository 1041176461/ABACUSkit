'''
Date: 2021-03-29 21:33:09
LastEditors: jiyuyang
LastEditTime: 2021-04-23 16:38:39
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.utils.tools import read_json
from pyautotest.utils.typings import *
from pyautotest.postprocess.plot import BandPlot
from pyautotest.calculations.structure import read_kpt

import os
from typing import Sequence, Union

class Show:
    """Show auto-test information"""

    @classmethod
    def show_libinfo(cls, src: str_PathLike):
        """Show example library information
    
        :params src: path of library
        """

        print("--------------------------Library Information--------------------------", flush=True)
        print(f"Path: {src}", flush=True)
        print("Directory Structure:", flush=True)
        for index, subsrc in enumerate(os.listdir(src)):
            if os.path.isdir(os.path.join(src, subsrc)) and subsrc != "OUT.test":
                line = f" ({index+1}) " + subsrc
                print(f"{line}", flush=True)
                subline = '\t' + ', '.join(os.listdir(os.path.join(src, subsrc)))
                print(subline, flush=True)
            elif os.path.isfile(os.path.join(src, subsrc)):
                subline = '\t' + ', '.join(os.listdir(src))
                print(subline, flush=True)
                break
            else:
                raise FileNotFoundError(f"No information to show!")

    @classmethod
    def show_bandinfo(cls, datafile:Union[str_PathLike, muti_Path], kptfile:str_PathLike, efermi:Union[float, Sequence[float]]=None, energy_range:Sequence[float]=[], blabel:Union[str, Sequence[str]]=None, color:Union[str, Sequence[str]]=None, outfile:str_PathLike="band.png"):
        """Show band structure information
        
        :params datafile: path of band date file 
        :params kptfile: path of k-points file
        :params efermi: Fermi levels in unit eV, its length equals to `filename`
        :params energy_range: range of energy to plot, its length equals to two
        :params blabel: band labels, its length equals to `filename`.
        :params color: band colors, its length equals to `filename`.
        :params outfile: band picture file name. Default: 'band.png'
        """

        if isinstance(datafile, (str, PathLike)):
            BandPlot.singleplot(datafile, kptfile, efermi, energy_range, blabel, color, outfile)
        elif isinstance(datafile, (list, tuple)):
            BandPlot.multiplot(datafile, kptfile, efermi, energy_range, blabel, color, outfile)

    @classmethod
    def show_cmdline(cls, args):
        if args.lib:
            text = read_json(args.lib)
            cls.show_libinfo(text["src"])
        
        if args.band:
            text = read_json(args.band)
            filename = text["filename"]
            kptfile = text["kptfile"]
            efermi = text.pop("efermi", 0.0)
            energy_range = text.pop("energy_range", []) 
            blabel = text.pop("blabel", None) 
            color = text.pop("color", None) 
            outfile = text.pop("outfile", "band.png")
            cls.show_bandinfo(filename, kptfile, efermi, energy_range, blabel, color, outfile)
