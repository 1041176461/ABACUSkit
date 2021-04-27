'''
Date: 2021-03-10 22:20:20
LastEditors: jiyuyang
LastEditTime: 2021-04-18 16:58:56
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.calculations.baseclass import ABACUSCalculation

import glob
import shutil
from pathlib import Path

class NSCF(ABACUSCalculation):
    """NSCF calculation"""

    def __init__(self, input_dict, src, **kwargs) -> None:
        """Set input parameters of nscf calcultion
        
        :params input_dict: dict of input parameters
        :params src: path of example which will be tested
        """
        
        super().__init__(input_dict, src, **kwargs)
        self.input_dict["calculation"] = "nscf"

class BAND(NSCF):
    """Band calculation"""

    def __init__(self, input_dict, src, **kwargs) -> None:
        """Set input parameters of band calcultion"""

        super().__init__(input_dict, src, **kwargs)
        self.input_dict["start_charge"] = "file"
        self.input_dict["out_band"] = 1

        if "density_file" in kwargs.keys():
            self.density_file = kwargs["density_file"]
        else:
            self.density_file = []

    def _prepare(self, dst, **kwargs):
        """Prepare input files for nscf calculation e.g. INPUT, STRU, KPT, orbital and pseudopotential files"""

        super()._prepare(dst, **kwargs)
        
        # KPT check
        with open(Path(dst, "KPT"), 'r') as file:
            if "Line\n" not in file.readlines():
                raise Exception(f"The `KPT` file in {dst} is not proper to band calculation")

        # Density file
        outdir = Path(dst, "OUT.test")
        if not self.density_file:
            if not glob.glob(str(outdir)+"/SPIN*_CHG") or not glob.glob(str(outdir)+"/HR_exx_*"):
                raise FileNotFoundError(f"Can not find density files in {outdir}")
        else:
            if not outdir.exists:
                outdir.mkdir()
            for i in self.density_file:
                shutil.copyfile(Path(self.src, i), Path(outdir, i))

    def _parse(self, dst, **kwargs):
        """parse output of nscf calculation
        
        :param dst: directory where calculation executes
        :return: if `BANDS_*.dat` exists, return {"band_file" : 1}. If not, raise FileNotFoundError
        """
        
        res = {}
        outdir = Path(dst, "OUT.test")
        band = glob.glob(str(outdir)+"/BANDS_*.dat")
        if band:
            res["band_file"] = 1    # if `band_file` exists, set 1 for autotest check
        else:
            raise FileNotFoundError(f"`BANDS_*.dat` is not found in {outdir}")
        
        return res