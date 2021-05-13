'''
Date: 2021-05-08 11:47:09
LastEditors: jiyuyang
LastEditTime: 2021-05-08 11:47:09
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.utils.tools import list_elem2str
from pyautotest.calculations.structure import read_kpt
from pyautotest.utils.typings import *

import re
import numpy as np
from collections import OrderedDict, namedtuple
from typing import Sequence, Tuple
from matplotlib import axes
import matplotlib.pyplot as plt

def set_range(energy:Sequence, energy_range:Sequence=[]) -> tuple:
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

def energy_minus_efermi(energy:Sequence , efermi:float) -> np.ndarray:
    """Return energy after subtracting the Fermi level

    :params efermi: Fermi level in unit eV
    """

    return np.array(energy)-efermi


class BandPlot:
    """Plot band structure"""

    @classmethod
    def read(cls, filename:str_PathLike) -> Tuple[np.ndarray, np.ndarray]:
        """Read band data file and return k-points and energy
        
        :params filename: string of band data file
        """

        data = np.loadtxt(filename)
        X, y= np.split(data, (1, ), axis=1)
        x = X.flatten()
        return x, y

    @classmethod
    def _set_figure(cls, ax:axes.Axes, index:dict, range:Sequence):
        """set figure and axes for plotting
        
        :params ax: matplotlib.axes.Axes object
        :params index: dict of label of points of x-axis and its index in data file. Range of x-axis based on index.value()
        :params range: range of y-axis
        """

        keys = []
        values = []
        for t in index:
            if isinstance(t, tuple):
                keys.append(t[0])
                values.append(t[1])
            elif isinstance(t, (int, float)):
                keys.append('')
                values.append(t)

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
        ax.axhline(0, linestyle="--", c='b', lw=1.0)
        handles, labels = ax.get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())

    @classmethod
    def plot(cls, x:Sequence, y:Sequence, index:Sequence, efermi:float=0, energy_range:Sequence[float]=[], label:str=None, color:str=None, outfile:str_PathLike='band.png'):
        """Plot band structure
        
        :params x, y: x-axis and y-axis coordinates
        :params index: special k-points label and its index in data file
        :params efermi: Fermi level in unit eV
        :params energy_range: range of energy to plot, its length equals to two
        :params label: band label. Default: ''
        :params color: band color. Default: 'black'
        :params outfile: band picture file name. Default: 'band.png'
        """

        fig, ax = plt.subplots()

        if not color:
            color = 'black'

        kpoints, energy = x, y
        energy = energy_minus_efermi(energy, efermi)
        energy_range = set_range(energy, energy_range)

        ax.plot(kpoints, energy, lw=0.8, color=color, label=label)
        cls._set_figure(ax, index, energy_range)

        plt.savefig(outfile)

    @classmethod
    def singleplot(cls, datafile:str_PathLike, kptfile:str=[], efermi:float=0, energy_range:Sequence[float]=[], label:str=None, color:str=None, outfile:str_PathLike='band.png'):
        """Plot band structure using data file
        
        :params datafile: string of band date file
        :params kptfile: k-point file
        :params efermi: Fermi level in unit eV
        :params energy_range: range of energy to plot, its length equals to two
        :params label: band label. Default: ''
        :params color: band color. Default: 'black'
        :params outfile: band picture file name. Default: 'band.png'
        """

        fig, ax = plt.subplots()
        kpt = read_kpt(kptfile)

        if not color:
            color = 'black'

        kpoints, energy = cls.read(datafile)
        energy = energy_minus_efermi(energy, efermi)
        energy_range = set_range(energy, energy_range)

        ax.plot(kpoints, energy, lw=0.8, color=color, label=label)
        cls.info(kpt.kpath, energy)
        index = kpt.label_special_k
        cls._set_figure(ax, index, energy_range)

        plt.savefig(outfile)
        
    @classmethod
    def multiplot(cls, datafile:muti_Path, kptfile:str=[], efermi:Sequence[float]=[], energy_range:Sequence[float]=[], label:Sequence[str]=None, color:Sequence[str]=None, outfile:str_PathLike='band.png'):
        """Plot more than two band structures using data file
        
        :params datafile: list of path of band date file 
        :params kptfile: k-point file
        :params efermi: list of Fermi levels in unit eV, its length equals to `filename`
        :params energy_range: range of energy to plot, its length equals to two
        :params label: list of band labels, its length equals to `filename`.
        :params color: list of band colors, its length equals to `filename`.
        :params outfile: band picture file name. Default: 'band.png'
        """

        fig, ax = plt.subplots()
        kpt= read_kpt(kptfile)

        if not efermi:
            efermi = [0.0 for i in range(len(datafile))]
        if not label:
            label = ['' for i in range(len(datafile))]
        if not color:
            color = ['black' for i in range(len(datafile))]

        emin = -np.inf
        emax = np.inf
        for i, file in enumerate(datafile):
            kpoints, energy = cls.read(file)
            energy = energy_minus_efermi(energy, efermi[i])
            energy_min = np.min(energy)
            energy_max = np.max(energy)
            if energy_min>emin:
                emin = energy_min
            if energy_max<emax:
                emax = energy_max

            ax.plot(kpoints, energy, lw=0.8, color=color[i], label=label[i])
            cls.info(kpt.kpath, energy)
        
        energy_range = set_range([emin, emax], energy_range)

        index = kpt.label_special_k
        cls._set_figure(ax, index, energy_range)

        plt.savefig(outfile)

    @classmethod
    def bandgap(cls, energy:Sequence):
        """Calculate band gap
        
        :params energy: band energy after subtracting the Fermi level
        """

        e_T = energy.T
        num_gt_Ef = (e_T > 0).sum(axis=1)

        Band = namedtuple('Band', ['band_index', 'band', 'value', 'k_index'])

        # valance band
        band_vbm_index = np.where(num_gt_Ef == 0)[0]
        band_vbm = e_T[band_vbm_index]
        evbm = np.max(band_vbm)
        k_vbm_index = np.where(band_vbm == evbm)[1]
        vbm = Band(band_vbm_index, band_vbm, evbm, k_vbm_index)

        # conduct band
        band_cbm_index = np.where(num_gt_Ef != 0)[0]
        band_cbm= e_T[band_cbm_index]
        ecbm = np.min(band_cbm)
        k_cbm_index = np.where(band_cbm == ecbm)[1]
        cbm = Band(band_cbm_index, band_cbm, ecbm, k_cbm_index)

        gap = ecbm-evbm

        return gap, vbm, cbm

    @classmethod
    def info(cls, kpath:Sequence, energy:Sequence):
        """Output the information of band structure
        
        :params kpath: k-points path 
        :params energy: band energy after subtracting the Fermi level
        """

        def band_type(vbm_x, cbm_x):
            longone, shortone = (vbm_x, cbm_x) if len(vbm_x) >= len(cbm_x) else (cbm_x, vbm_x)
            for i in shortone:
                if i in longone:
                    btype = "Direct"
                else:
                    btype = "Indirect"
            return btype

        gap, vbm, cbm = cls.bandgap(energy)
        print("--------------------------Band Structure--------------------------", flush=True)
        print(f"{'Band character:'.ljust(30)}{band_type(vbm.k_index, cbm.k_index)}", flush=True)
        print(f"{'Band gap(eV):'.ljust(30)}{gap: .4f}", flush=True)
        print(f"{'Band index:'.ljust(30)}{'HOMO'.ljust(10)}{'LUMO'}", flush=True)
        print(f"{''.ljust(30)}{str(vbm.band_index[-1]).ljust(10)}{str(cbm.band_index[0])}", flush=True)
        print(f"{'Eigenvalue of VBM(eV):'.ljust(30)}{vbm.value: .4f}", flush=True)
        print(f"{'Eigenvalue of CBM(eV):'.ljust(30)}{cbm.value: .4f}", flush=True)
        vbm_k = np.unique(kpath[vbm.k_index], axis=0)
        cbm_k = np.unique(kpath[cbm.k_index], axis=0)
        print(f"{'Location of VBM'.ljust(30)}{' '.join(list_elem2str(vbm_k[0]))}", flush=True)
        for i, j in enumerate(vbm_k):
            if i != 0:
                print(f"{''.ljust(30)}{' '.join(list_elem2str(j))}", flush=True)
        print(f"{'Location of CBM'.ljust(30)}{' '.join(list_elem2str(cbm_k[0]))}", flush=True)
        for i, j in enumerate(cbm_k):
            if i != 0:
                print(f"{''.ljust(30)}{' '.join(list_elem2str(j))}", flush=True)

class DosPlot:
    """Plot density of state(DOS)"""

    @classmethod
    def _val_to_zero(cls, string):
        if eval(string) <= 1e-20:
            string = "0"
        return string

    @classmethod
    def read(cls, tdosfile:str_PathLike='', pdosfile:str_PathLike='') -> tuple:
        """Read DOS data file, if both `tdosfile` and `pdosfile` set, it will ony read `tdosfile`
        
        :params tdosfile: string of TDOS data file
        :params pdosfile: string of PDOS data file
        """

        if tdosfile:
            dosdata = np.loadtxt(tdosfile)
            nsplit = dosdata.shape[1]
            return np.split(dosdata, nsplit, axis=1)
        
        elif pdosfile:
            e_list = []
            orbitals = []
            with open(pdosfile, "r") as f:
                for line in f:
                    if re.match(r"</pdos>", line):
                        break
                    elif re.match(r"<nspin>[0-9]</nspin>", line):
                        nspin = int(re.match(r"(<nspin>)([0-9])(</nspin>)", line).group(2))
                        continue
                    elif re.match(r"<norbitals>26</norbitals>", line):
                        norbitals = int(re.match(r"(<norbitals>)([0-9]+)(</norbitals>)", line).group(2))
                        continue
                    elif re.match(r"<energy_values units=\"[a-zA-Z]+\">", line):
                        e_unit = re.match(r"(<energy_values units=\")([a-zA-Z]+)(\">)", line).group(2)
                        for line in f:
                            if re.match("</energy_values>", line):
                                break
                            e = float(re.compile(r"\s+").split(line, maxsplit=1)[1].strip())
                            e_list.append(e)
                        npoints = len(e_list)
                        continue
                    elif re.match(r"<orbital", line):
                        orb = OrderedDict()
                        for line in f:
                            if re.match(r"\s*index=\"\s*[0-9]+\"", line):
                                orb['index'] = int(re.compile(r"\s+").split(line, maxsplit=2)[-1].strip("\"\n"))
                                continue
                            elif re.match(r"\s*atom_index=\"\s*[0-9]+\"", line):
                                orb['atom_index'] = int(re.compile(r"\s+").split(line, maxsplit=2)[-1].strip("\"\n"))
                                continue
                            elif re.match(r"\s*species=\"\s*[A-Za-z]+\"", line):
                                orb['species'] = re.compile(r"\=").split(line, maxsplit=1)[-1].strip("\"\n")
                                continue
                            elif re.match(r"\s*l=\"\s*[0-9]+\"", line):
                                orb['l'] = int(re.compile(r"\s+").split(line, maxsplit=2)[-1].strip("\"\n"))
                                continue
                            elif re.match(r"\s*m=\"\s*[0-9]+\"", line):
                                orb['m'] = int(re.compile(r"\s+").split(line, maxsplit=2)[-1].strip("\"\n"))
                                continue
                            elif re.match(r"\s*z=\"\s*[0-9]+\"", line):
                                orb['z'] = int(re.compile(r"\s+").split(line, maxsplit=2)[-1].strip("\"\n"))
                                continue
                            elif re.match(r"<data>", line):
                                orb['data'] = np.zeros((npoints, nspin), dtype=np.float32)
                                for i, line in enumerate(f):
                                    if re.match(r"</data>", line):
                                        orbitals.append(orb)
                                        orb = OrderedDict()
                                        break
                                    val = re.compile(r"\s").split(line)
                                    while '' in val:
                                        val.remove('')
                                    orb['data'][i] = [j for j in map(cls._val_to_zero, val)]
                                continue
                        continue
            return np.reshape(e_list, newshape=(-1, 1)), orbitals

    @classmethod
    def _set_figure(cls, ax:axes.Axes, dos_range:Sequence, energy_range:Sequence):
        """set figure and axes for plotting
        
        :params ax: matplotlib.axes.Axes object
        :params dos_range: range of dos
        :params energy_range: range of energy
        """

        # x-axis
        ax.set_xlim(dos_range[0], dos_range[1])
        ax.set_xlabel("DOS")

        # y-axis
        ax.set_ylim(energy_range[0], energy_range[1])
        ax.set_ylabel(r"$E-E_{fermi}(eV)$")

        # others
        ax.axhline(0, linestyle="--", c='b', lw=1.0)
        ax.legend()
    
    @classmethod
    def _plot(cls, res:tuple, efermi:float=0, energy_range:Sequence[float]=[]):

        fig, ax = plt.subplots()

        nsplit = len(res)
        if nsplit == 2:
            energy, dos = res
            energy = energy_minus_efermi(energy, efermi)
            dos_range = (np.min(dos), np.max(dos))
            energy_range = set_range(energy, energy_range)
            ax.plot(dos, energy, lw=0.8, c='gray', linestyle='-', label='TDOS')
                
        elif nsplit == 3:
            energy, dos_up, dos_dw = res
            energy = energy_minus_efermi(energy, efermi)
            dos_dw = -dos_dw
            dos_range = (np.min(dos_dw), np.max(dos_up))
            energy_range = set_range(energy, energy_range)
            ax.plot(dos_up, energy, lw=0.8, c='gray', linestyle='-', label=r'$TDOS \uparrow$')
            ax.plot(dos_dw, energy, lw=0.8, c='gray', linestyle='--', label=r'$TDOS \downarrow$')
        
        cls._set_figure(ax, dos_range, energy_range)

        return ax

    @classmethod
    def _where_sum(cls, orbitals, key, value):
        nsplit = orbitals[0]["data"].shape[1]
        res = np.zeros_like(orbitals[0]["data"], dtype=np.float32)
        for orb in orbitals:
            if orb[key] == value:
                res = res + orb['data']
        return res, nsplit

    @classmethod
    def _all_sum(cls, orbitals):
        nsplit = orbitals[0]["data"].shape[1]
        res = np.zeros_like(orbitals[0]["data"], dtype=np.float32)
        for orb in orbitals:
            res = res + orb['data']
        return res, nsplit

    @classmethod
    def plot(cls, tdosfile:str_PathLike='', pdosfile:str_PathLike='', efermi:float=0, energy_range:Sequence[float]=[], choose:dict={}, outfile:str_PathLike='dos.png'):
        """Plot total dos or partial dos, if both `tdosfile` and `pdosfile` set, it will ony read `tdosfile`
        
        :params tdosfile: string of TDOS data file
        :params pdosfile: string of PDOS data file
        :params efermi: Fermi level in unit eV
        :params energy_range: range of energy to plot, its length equals to two
        :params choose: specify label for plotting PDOS. Supported keys: 'atom_index', 'species', 'l', 'm', 'z'
        """

        res = cls.read(tdosfile, pdosfile)

        if tdosfile:
            cls._plot(res, efermi, energy_range)
        elif pdosfile:
            energy, orbitals = res
            # TDOS
            dos, nsplit = cls._all_sum(orbitals)
            if nsplit == 1:
                ax = cls._plot((energy, dos), efermi, energy_range)
            elif nsplit == 2:
                dos_up, dos_dw = np.split(dos, nsplit, axis=1)
                ax = cls._plot((energy, dos_up, dos_dw), efermi, energy_range)
            # PDOS
        #    if choose:
        #        xmax, xmin = [], []
        #        energy = energy_minus_efermi(energy, efermi)
        #        for key, value in choose.items():
        #            for val in value:
        #                dos, nsplit = cls._where_sum(orbitals, key, val)
        #                if nsplit == 1:
        #                    ax.plot(dos, energy, lw=0.8, linestyle='-', label=f'{val}')
        #                elif nsplit == 2:
        #                    dos_up, dos_dw = np.split(dos, nsplit, axis=1)
        #                    ax.plot(dos_up, energy, lw=0.8, linestyle="-", label=f"{val}"+r"$\uparrow$")
        #                    ax.plot(dos_dw, energy, lw=0.8, linestyle="--", label=f"{val}"+r"$\downarrow$")
        #            xmax.append(np.max(dos))
        #            xmin.append(np.min(dos))
        #        dos_range = (np.min(xmin), np.max(xmax))
        #        energy_range = set_range(energy, energy_range)
        #        cls._set_figure(ax, dos_range, energy_range)

        plt.savefig(outfile)