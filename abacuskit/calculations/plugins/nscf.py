'''
Date: 2021-03-10 22:20:20
LastEditors: jiyuyang
LastEditTime: 2021-04-18 16:58:56
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from abacuskit.calculations.baseclass import ABACUSCalculation
from abacuskit.calculations.structure import Stru, Kpt
from abacuskit.utils.typings import *

import glob
import shutil
from pathlib import Path

class NSCF(ABACUSCalculation):
    """NSCF calculation"""

    def __init__(self, input_dict: dict, stru: typing.Optional[Stru], kpt: typing.Optional[Kpt], **kwargs) -> None:
        """Set input parameters of nscf calcultion
        
        :params input_dict: dict of input parameters
        :params stru: object of `abacuskit.calculations.structure.Stru`
        :params kpt: object of `abacuskit.calculations.structure.Kpt`
        """

        super().__init__(input_dict, stru, kpt, **kwargs)
        self.input_dict["calculation"] = "nscf"

class BAND(NSCF):
    """Band calculation"""

    def __init__(self, input_dict: dict, stru: typing.Optional[Stru], kpt: typing.Optional[Kpt], density_file: str_PathLike="", **kwargs) -> None:
        """Set input parameters of band calcultion
        
        :params input_dict: dict of input parameters
        :params stru: object of `abacuskit.calculations.structure.Stru`
        :params kpt: object of `abacuskit.calculations.structure.Kpt`
        """

        super().__init__(input_dict, stru, kpt, **kwargs)
        self.input_dict["start_charge"] = "file"
        self.input_dict["out_band"] = 1
        self.density_file = density_file  # should be absolute path of density file

    def _prepare(self, **kwargs):
        """Prepare input files for nscf calculation"""
        
        super()._prepare(**kwargs)
        
        # KPT check
        if self.kpt.mode != "Line":
            raise Exception(f"The `KPT` file is not proper to band calculation")

        # Density
        outdir = "OUT.test"
        Path(outdir).mkdir(exist_ok=True)
        if not glob.glob(outdir+"/SPIN*_CHG") or not glob.glob(outdir+"/HR_exx_*"):
            for i in self.density_file:
                shutil.copy(i, outdir)

    def _parse(self, **kwargs) -> dict:
        """parse output of nscf calculation
        
        :return: if `BANDS_*.dat` exists, return {"band_file" : 1}. If not, raise FileNotFoundError
        """
        
        res = {}
        outdir = "OUT.test"
        band = glob.glob(outdir+"/BANDS_*.dat")
        if band:
            res["band_file"] = 1    # if `band_file` exists, set 1 for autotest check
        else:
            raise FileNotFoundError(f"`BANDS_*.dat` is not found in {outdir}")
        
        return res