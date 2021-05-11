'''
Date: 2021-05-08 11:47:09
LastEditors: jiyuyang
LastEditTime: 2021-05-08 11:47:09
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.utils.typings import *

import numpy as np
from collections import OrderedDict
from typing import Sequence, Tuple
from matplotlib import axes
import matplotlib.pyplot as plt

class BandPlot:
    """Plot band structure"""

    def __init__(self) -> None:
        """Set figure"""

    @classmethod
    def read(cls, filename:str) -> Tuple[np.ndarray, np.ndarray]:
        """Read band date file and return k-points and energy
        
        :params filename: string of band date file
        """

        data = np.loadtxt(filename)
        X, y= np.split(data, (1, ), axis=1)
        x = X.flatten()
        return x, y
    
    @classmethod
    def _set_range(cls, energy:Sequence, energy_range:Sequence=[]) -> tuple:
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
    def energy_minus_efermi(cls, energy:Sequence , efermi:float) -> np.ndarray:
        """Return energy after subtracting the Fermi level
        
        :params efermi: Fermi level in unit eV
        """

        return np.array(energy)-efermi

    @classmethod
    def _set_figure(cls, ax:axes.Axes, index:dict, range:Sequence):
        """set figure and axes for plotting
        
        :params ax: matplotlib.axes.Axes object
        :params index: dict of label of points of x-axis and its index in data file. Range of x-axis based on index.value()
        :params range: range of y-axis
        """

        keys = []
        values = []
        for key, value in index:
            keys.append(key)
            values.append(value)

        # x-axis
        ax.set_xticks(values)
        ax.set_xticklabels(keys)
        ax.set_xlim(values[0], values[-1])
        ax.set_xlabel("Wave Vector")

        # y-axis
        ax.set_ylim(range[0], range[1])
        ax.set_ylabel(r"$E-E_{fermi}(eV)$")

        # others
        ax.grid(axis='x', lw=1.2)
        ax.axhline(0, linestyle="--", lw=1.0)
        handles, labels = ax.get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())

    @classmethod
    def plot(cls, filename:str_PathLike, index:dict, efermi:float=0, energy_range:Sequence[float]=[], label:str=None, color:str=None, outfile:str_PathLike='band.png'):
        """Plot band structure
        
        :params filename: string of band date file
        :params index: dict of special k-points label and its index in data file
        :params efermi: Fermi level in unit eV
        :params energy_range: range of energy to plot, its length equals to two
        :params label: band label. Default: ''
        :params color: band color. Default: 'black'
        :params outfile: band picture file name. Default: 'band.png'
        """

        fig, ax = plt.subplots()

        if not color:
            color = 'black'

        kpoints, energy = cls.read(filename)
        energy = cls.energy_minus_efermi(energy, efermi)
        energy_range = cls._set_range(energy, energy_range)

        ax.plot(kpoints, energy, lw=0.8, color=color, label=label)
        cls._set_figure(ax, index, energy_range)

        plt.savefig(outfile)
        
    @classmethod
    def multiplot(cls, filename:muti_Path, index:dict, efermi:Sequence[float]=[], energy_range:Sequence[float]=[], label:Sequence[str]=None, color:Sequence[str]=None, outfile:str_PathLike='band.png'):
        """Plot more than two band structures
        
        :params filename: list of path of band date file 
        :params index: list of special k-points label and its index in data file e.g. [("G", 1), ("L", 10), ...]
        :params efermi: list of Fermi levels in unit eV, its length equals to `filename`
        :params energy_range: range of energy to plot, its length equals to two
        :params label: list of band labels, its length equals to `filename`.
        :params color: list of band colors, its length equals to `filename`.
        :params outfile: band picture file name. Default: 'band.png'
        """
        
        fig, ax = plt.subplots()

        if not efermi:
            efermi = [0.0 for i in range(len(filename))]
        if not label:
            label = ['' for i in range(len(filename))]
        if not color:
            color = ['black' for i in range(len(filename))]

        emin = -np.inf
        emax = np.inf
        for i, file in enumerate(filename):
            kpoints, energy = cls.read(file)
            energy = cls.energy_minus_efermi(energy, efermi[i])
            energy_min = np.min(energy)
            energy_max = np.max(energy)
            if energy_min>emin:
                emin = energy_min
            if energy_max<emax:
                emax = energy_max

            ax.plot(kpoints, energy, lw=0.8, color=color[i], label=label[i])
        
        energy_range = cls._set_range([emin, emax], energy_range)
        cls._set_figure(ax, index, energy_range)

        plt.savefig(outfile)