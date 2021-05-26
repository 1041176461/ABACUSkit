'''
Date: 2021-05-25 23:01:57
LastEditors: jiyuyang
LastEditTime: 2021-05-25 23:01:58
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from abacuskit.postprocess.symmetry import Spacegroup
from abacuskit.utils.IO import read_stru
from pathlib import Path

ntype = 4
filename = r"C:\Users\jiyuyang\Desktop\STRU"
pardir = Path(filename).parent
stru = read_stru(4, r"C:\Users\jiyuyang\Desktop\STRU")
obj = Spacegroup(stru)
pri_stru, kpt = obj.get_kpath()
pri_stru.write_stru(pardir/"STRU-pri")
kpt.write_kpt(pardir/"KPT-band")