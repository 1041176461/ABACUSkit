'''
Date: 2021-03-29 21:35:30
LastEditors: jiyuyang
LastEditTime: 2021-04-23 16:51:39
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

import re
import sys
import json
import numpy as np
from collections import defaultdict

def read_json(filename=""):
    """ Read json file and return dict
    
    :params filename: json file
    """
    with open(filename, 'r') as file:
        text = json.load(file)
    return text

def write_json(filename="", new_filename="", **kwargs):
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

def read_cif(filename=""):
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

def add_path(path=""):
    """ Add a path into sys.path
    
    :params path: path which will be added
    """
    if path not in sys.path:
        sys.path.append(path)

def skip_notes(line=""):
    """Delete comments lines with '#' or '//'

    :params line: line will be handled 
    """
    line = re.compile(r"#.*").sub("",line)
    line = re.compile(r"//.*").sub("",line)
    line = line.strip()
    return line

def ignore_lines(file, n=0):
    """Ignore n lines in file
    
    :params file: file descriptor
    :params n: number of lines will be ignored
    """
    for i in range(n):
        file.readline()

def list_elem2strip(a=[]):
    """Strip element of list with `str` type"""
    def list_strip(s):
        return s.strip()
    return list(map(list_strip, a))

def search_sentence(file, sentence):
    """Search sentence in file
    
    :params file: file descriptor
    :params sentence: sentence will be searched
    """
    
    if isinstance(sentence, str):
        sentence = sentence.strip()
        for line in file:
            if skip_notes(line).strip() == sentence:
                return True
        return False
    elif isinstance(sentence, list):
        sentence = list_elem2strip(sentence)
        for line in file:
            if skip_notes(line).strip() in sentence:
                return True
        return False

def list_elem2str(a=[]):
    """Convert type of list element to str
    
    :params a: 1-D list
    """
    return list(map(str, a))

def folder_name(T1, T2, i_dis):
    return f"{T1}-{T2}_{i_dis}"