'''
Date: 2021-05-14 16:03:47
LastEditors: jiyuyang
LastEditTime: 2021-05-14 16:10:36
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''
from pyautotest.calculations.structure import read_stru
from pathlib import Path
from pyautotest.utils.constants import BOHR_TO_A

tpath = "./" # target directory
ntype = 4

E = []
a = []
V = []

lat = 10.4805927277

for i in range(0, 5):
    spath = Path(tpath, f"{i}")
    obj = read_stru(ntype, Path(spath, "STRU"))
    obj.set_energy(Path(spath, "OUT.test/running_scf.log"))
    V.append(obj.volume)
    E.append(obj.energy)
    a.append(obj.lat0*lat*BOHR_TO_A)

with open("V_E_res", 'w') as f1, open("a_E_res", 'w') as f2:
    for i in range(0, 5):
        f1.write(f"{V[i]}\t{E[i]}\n")
        f2.write(f"{a[i]}\t{E[i]}\n")

import os
factor = 1
fit_script = "fit-bulk.py"
os.system(f"python {fit_script} V_E_res > hse_V_fit")
os.system(f"python {fit_script} -l {factor} a_E_res > a_fit")
os.system(f"python {fit_script} -p V_E_res")