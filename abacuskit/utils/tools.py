'''
Date: 2021-03-29 21:35:30
LastEditors: jiyuyang
LastEditTime: 2021-05-16 14:30:15
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

import re
import string
import sys
from typing import List, Union

from abacuskit.utils.typings import *


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
    line = re.compile(r"#.*").sub("", line)
    line = re.compile(r"//.*").sub("", line)
    line = line.strip()
    return line


def ignore_lines(file: str_PathLike, n: int = 0):
    """Ignore n lines in file

    :params file: file descriptor
    :params n: number of lines will be ignored
    """
    for i in range(n):
        file.readline()


def list_elem2strip(a: List[str], ds=string.whitespace) -> List[str]:
    """Strip element of list with `str` type"""
    def list_strip(s):
        return s.strip(ds)
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


def list_elem_2float(a: List[str]) -> List[float]:
    """Convert type of list element to float

    :params a: 1-D list
    """
    return list(map(float, a))


def list_elem_2int(a: List[str]) -> List[int]:
    """Convert type of list element to int

    :params a: 1-D list
    """
    return list(map(int, a))


def folder_name(T1: str, T2: str, i_dis: Union[int, float, str]) -> str:
    return f"{T1}-{T2}_{i_dis}"


def delete_key(input_dict):
    key_list = ["ocp", "ocp_set", "nelec", "out_charge"]
    for i in key_list:
        input_dict.pop(i, None)


def get_input_line(input_dict: dict):
    """Return input lines in INPUT file

    :params input_dict: dict of input parameters
    """

    lines = []
    lines.append("INPUT_PARAMETERS")
    for key, value in input_dict.items():
        if value:
            lines.append(f"{key.ljust(30)}{value}")
    return '\n'.join(lines)
