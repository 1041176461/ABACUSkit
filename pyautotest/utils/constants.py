'''
Date: 2021-05-12 00:33:51
LastEditors: jiyuyang
LastEditTime: 2021-05-12 00:33:51
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from scipy.constants import physical_constants

BOHR_TO_A = physical_constants["atomic unit of length"][0]/physical_constants["Angstrom star"][0]
Hartree_TO_eV = physical_constants["Hartree energy in eV"][0]
Rydberg_TO_eV = physical_constants["Rydberg constant times hc in eV"]