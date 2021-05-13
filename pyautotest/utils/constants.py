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

def get_angular_momentum_label(index:int) -> str: 
    """Atomic orbital angular momentum label from index
    
    :params index: 0 or 1 or 2 or 3
    """

    if index == 0:
        return 's'
    elif index == 1:
        return 'p'
    elif index == 2:
        return 'd'
    elif index == 3:
        return 'f'

def get_angular_momentum_index(label:str) -> int: 
    """Atomic orbital angular momentum index from label
    
    :params label: 's' or 'p' or 'd' or 'f'
    """

    if label == 's':
        return 0
    elif label == 'p':
        return 1
    elif label == 'd':
        return 2
    elif label == 'f':
        return 3