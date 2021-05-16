'''
Date: 2021-05-08 11:47:09
LastEditors: jiyuyang
LastEditTime: 2021-05-08 11:47:09
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.utils.constants import get_angular_momentum_label
from pyautotest.utils.tools import list_elem2str
from pyautotest.utils.IO import read_kpt
from pyautotest.utils.typings import *

import re
import numpy as np
from collections import OrderedDict, namedtuple
from typing import Dict, Sequence, Tuple, Union, List
from matplotlib import axes
import matplotlib.pyplot as plt
from functools import lru_cache

def energy_minus_efermi(energy:Sequence , efermi:float) -> np.ndarray:
    """Return energy after subtracting the Fermi level

    :params efermi: Fermi level in unit eV
    """

    return np.array(energy)-efermi


class BandPlot:
    """Plot band structure"""

    @classmethod
    @lru_cache(maxsize=None, typed=False)
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
        if range:
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
    @lru_cache(maxsize=None, typed=False)
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
    def _set_figure(cls, ax:axes.Axes, energy_range:Sequence, dos_range:Sequence):
        """set figure and axes for plotting
        
        :params ax: matplotlib.axes.Axes object
        :params dos_range: range of dos
        :params energy_range: range of energy
        """

        # y-axis
        if dos_range:
            ax.set_ylim(dos_range[0], dos_range[1])
        ax.set_ylabel("DOS")

        # x-axis
        if energy_range:
            ax.set_xlim(energy_range[0], energy_range[1])
        ax.set_xlabel(r"$E-E_{fermi}(eV)$")

        # others
        ax.axvline(0, linestyle="--", c='b', lw=1.0)
        ax.legend()
    
    @classmethod
    def _tplot(cls, res:tuple, efermi:float=0, energy_range:Sequence[float]=[], dos_range:Sequence[float]=[]):

        fig, ax = plt.subplots()

        nsplit = len(res)
        if nsplit == 2:
            energy, dos = res
            energy = energy_minus_efermi(energy, efermi)
            ax.plot(energy, dos, lw=0.8, c='gray', linestyle='-', label='TDOS')
                
        elif nsplit == 3:
            energy, dos_up, dos_dw = res
            energy = energy_minus_efermi(energy, efermi)
            dos_dw = -dos_dw
            ax.plot(energy, dos_up, lw=0.8, c='gray', linestyle='-', label=r'$TDOS \uparrow$')
            ax.plot(energy, dos_up, lw=0.8, c='gray', linestyle='--', label=r'$TDOS \downarrow$')

        return ax, energy_range, dos_range

    @classmethod
    def _all_sum(cls, orbitals:dict) -> Tuple[np.ndarray, int]:
        nsplit = orbitals[0]["data"].shape[1]
        res = np.zeros_like(orbitals[0]["data"], dtype=np.float32)
        for orb in orbitals:
            res = res + orb['data']
        return res, nsplit

    @classmethod
    def plot(cls, tdosfile:str_PathLike='', pdosfile:str_PathLike='', efermi:float=0, energy_range:Sequence[float]=[], dos_range:Sequence[float]=[], species:Union[Sequence[str], Dict[str, List[int]]]=[], tdosfig:str_PathLike='tdos.png', pdosfig:str_PathLike='pdos.png'):
        """Plot total dos or partial dos, if both `tdosfile` and `pdosfile` set, it will ony read `tdosfile`
        
        :params tdosfile: string of TDOS data file
        :params pdosfile: string of PDOS data file
        :params efermi: Fermi level in unit eV
        :params energy_range: range of energy to plot, its length equals to two
        :params dos_range: range of dos to plot, its length equals to two
        :params species: list of atomic species or dict of atomic species and its angular momentum list
        """

        res = cls.read(tdosfile, pdosfile)

        if tdosfile:
            ax, energy_range, dos_range = cls._tplot(res, efermi, energy_range, dos_range)
            cls._set_figure(ax, energy_range, dos_range)
            plt.savefig(tdosfig)

        elif pdosfile and species:
            elements = []
            momentum = []
            if isinstance(species, (list, tuple)):
                elements = species
            elif isinstance(species, dict):
                elements = list(species.keys())
                momentum = list(species.values())
            if not elements:
                raise TypeError("Only when `pdosfile` and `species` are both set, it will plot PDOS.")
            energy, orbitals = res
            energy_f = energy_minus_efermi(energy, efermi)

            # TDOS
            dos, nsplit = cls._all_sum(orbitals)
            if nsplit == 1:
                ax, energy_range, dos_range = cls._tplot((energy, dos), efermi, energy_range, dos_range)
            elif nsplit == 2:
                dos_up, dos_dw = np.split(dos, nsplit, axis=1)
                ax, energy_range, dos_range = cls._tplot((energy, dos_up, dos_dw), efermi, energy_range, dos_range)
            for elem in elements:
                dos = np.zeros_like(orbitals[0]["data"], dtype=np.float32)
                for orb in orbitals:
                    if orb["species"] == elem:
                        dos += orb["data"]
                if nsplit == 1:
                    ax.plot(energy_f, dos, lw=0.8, linestyle='-', label=f'{elem}')
                elif nsplit == 2:
                    dos_up, dos_dw = np.split(dos, nsplit, axis=1)
                    ax.plot(energy_f, dos_up, lw=0.8, linestyle="-", label=f"{elem}"+r"$\uparrow$")
                    dos_dw = -dos_dw
                    ax.plot(energy_f, dos_dw, lw=0.8, linestyle="--", label=f"{elem}"+r"$\downarrow$")
            
            cls._set_figure(ax, energy_range, dos_range)
            plt.savefig(tdosfig)
            
            # PDOS            
            if momentum:
                fig, ax = plt.subplots(len(elements), 1, sharex=True, sharey=True)
                plt.subplots_adjust(hspace=0)
                for i, elem in enumerate(elements):
                    for j in momentum[i]:
                        dos = np.zeros_like(orbitals[0]["data"], dtype=np.float32)
                        for orb in orbitals:
                            if orb["species"] == elem and orb["l"] == j:
                                dos += orb["data"]
                        if nsplit == 1:
                            ax[i].plot(energy_f, dos, lw=0.8, linestyle='-', label=f'{elem}-{get_angular_momentum_label(j)}')
                        elif nsplit == 2:
                            dos_up, dos_dw = np.split(dos, nsplit, axis=1)
                            ax[i].plot(energy_f, dos_up, lw=0.8, linestyle="-", label=f"{elem}-{get_angular_momentum_label(j)}"+r"$\uparrow$")
                            dos_dw = -dos_dw
                            ax[i].plot(energy_f, dos_dw, lw=0.8, linestyle="--", label=f"{elem}-{get_angular_momentum_label(j)}"+r"$\downarrow$")

                    cls._set_figure(ax[i], energy_range, dos_range)

                plt.savefig(pdosfig)