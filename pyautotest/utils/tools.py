'''
Date: 2021-03-29 21:35:30
LastEditors: jiyuyang
LastEditTime: 2021-04-29 13:38:49
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.utils.typings import *

import re
import sys
import json
from collections import defaultdict
from typing import Tuple, List, Union

def read_json(filename: str_PathLike) -> dict:
    """ Read json file and return dict
    
    :params filename: json file
    """
    with open(filename, 'r') as file:
        text = json.load(file)
    return text

def write_json(filename: str_PathLike, new_filename: str_PathLike, **kwargs) -> str_PathLike:
    """ Read json file and modify some key-values, then write it to a new json file

    :params filename: json file to be read
    :params new_filename: json file to be written
    :params **kwargs: any key-value to be written to new_filename
    """
    with open(filename, 'r') as file:
        text = json.load(file)
        for key, value in kwargs.items():
            text[key] = value
    with open(new_filename, 'w') as file:
        json.dump(text, file, indent=4)
    return new_filename

def read_cif(filename: str_PathLike) -> Tuple[tuple, dict]:
    """Read cif file, return lattice and position
    
    :params filename: cif file
    """
    res = {}
    with open(filename, 'r') as file:
        for line in file:
            if re.search("_cell_length_a", line):
                a = float(line.split()[1])
            if re.search("_cell_length_b", line):
                b = float(line.split()[1])
            if re.search("_cell_length_c", line):
                c = float(line.split()[1])
            if re.search("_cell_angle_alpha", line):
                alpha = float(line.split()[1])
            if re.search("_cell_angle_beta", line):
                beta = float(line.split()[1])
            if re.search("_cell_angle_gamma", line):
                gamma = float(line.split()[1])
            if re.search("_atom_site_fract_z", line):
                position = defaultdict(list)
                for atom in file:
                    elem, x, y, z = atom.split()
                    position[elem].append((float(x), float(y), float(z)))
        lattice = (a, b, c, alpha, beta, gamma)
    
    return lattice, position

def add_path(path: str_PathLike):
    """ Add a path into sys.path
    
    :params path: path which will be added
    """
    if path not in sys.path:
        sys.path.append(path)

def skip_notes(line: str) -> str:
    """Delete comments lines with '#' or '//'

    :params line: line will be handled 
    """
    line = re.compile(r"#.*").sub("",line)
    line = re.compile(r"//.*").sub("",line)
    line = line.strip()
    return line

def ignore_lines(file: str_PathLike, n: int=0):
    """Ignore n lines in file
    
    :params file: file descriptor
    :params n: number of lines will be ignored
    """
    for i in range(n):
        file.readline()

def list_elem2strip(a: List[str]) -> List[str]:
    """Strip element of list with `str` type"""
    def list_strip(s):
        return s.strip()
    return list(map(list_strip, a))

def search_sentence(file: str_PathLike, sentence: str) -> Union[bool, str]:
    """Search sentence in file
    
    :params file: file descriptor
    :params sentence: sentence will be searched
    """
    
    if isinstance(sentence, str):
        sentence = sentence.strip()
        for line in file:
            line = skip_notes(line).strip()
            if line == sentence:
                return line
    elif isinstance(sentence, list):
        sentence = list_elem2strip(sentence)
        for line in file:
            line = skip_notes(line).strip()
            if line in sentence:
                return line

    file.seek(0, 0)
    return False

def list_elem2str(a: Union[List[float], List[int]]) -> List[str]:
    """Convert type of list element to str
    
    :params a: 1-D list
    """
    return list(map(str, a))

def folder_name(T1: str, T2: str, i_dis: Union[int, float, str]) -> str:
    return f"{T1}-{T2}_{i_dis}"

def get_input_line(input_dict:dict):
    """Return input lines in INPUT file
        
    :params input_dict: dict of input parameters
    """

    lines = []
    lines.append("INPUT_PARAMETERS")
    for key, value in input_dict.items():
        if value:
            lines.append(f"{key.ljust(30)}{value}")
    return '\n'.join(lines)

def write_input(input_dict:dict, filename:str_PathLike="INPUT"):
    """Write INPUT file based on input_dict
        
    :params input_dict: dict of input parameters
    """

    with open(filename, 'w') as file:
        file.write(get_input_line(input_dict))