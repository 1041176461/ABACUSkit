'''
Date: 2021-05-08 11:47:09
LastEditors: jiyuyang
LastEditTime: 2021-05-08 11:47:09
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

import numpy as np
import pandas as pd
import plotly.express as px
from typing import Sequence

class BandPlot:
    """Plot band structure"""

    def __init__(self) -> None:
        """Set figure"""

    @classmethod
    def read(cls, filename:str):
        """Read band date file and return k-points and energy
        
        :params filename: string of band date file
        """

        data = np.loadtxt(filename)
        X, y= np.split(data, (1, ), axis=1)
        x = X.flatten()
        return x, y
    
    @classmethod
    def set_range(cls, energy:Sequence, energy_range:Sequence=[]):
        """Set energy range
        
        :params energy: band energy list
        :params energy_range: range of energy to plot
        """

        length = len(energy_range)
        if length == 0:
            energy_range = (np.min(energy), np.max(energy))
        elif length == 1:
            energy_range = (np.min(energy), energy_range[0])
        elif length == 2:
            energy_range = (energy_range)
        else:
            raise ValueError("Length of `energy_range` must be less than or equal to 2.")

        return energy_range

    @classmethod
    def energy_minus_efermi(cls, energy:Sequence , efermi:float):
        """Return energy after subtracting the Fermi level
        
        :params efermi: Fermi level in unit eV
        """

        return np.array(energy)-efermi

    @classmethod
    def plot(cls, filename:str, efermi:float, index:dict, energy_range:Sequence=[]):
        """Plot band structure
        
        :params filename: string of band date file
        :params efermi: Fermi level in unit eV
        :params index: dict of special k-points and its index
        :params energy_range: range of energy to plot
        """

        kpoints, energy = cls.read(filename)
        energy = cls.energy_minus_efermi(energy, efermi)
        energy_range = cls.set_range(energy, energy_range)
        

