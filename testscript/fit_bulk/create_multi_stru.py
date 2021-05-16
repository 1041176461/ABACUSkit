from pyautotest.calculations import structure
from pathlib import Path
import shutil
import os

ppath = "/home/jiyy/ABACUS/AutoTest/new/testlib/MgO" # source file directory
stru_path = Path(ppath, "STRU")
filelist = os.listdir(ppath)
filelist.remove("STRU")
tpath = "./" # target directory
ntype = 2

obj = structure.read_stru(ntype, stru_path)
latc = obj.lat0
delta = 0.01
for i in range(-1, 4):
    obj.lat0 = latc+i*delta
    spath = Path(tpath, str(i+1))
    spath.mkdir()
    obj.write_stru(Path(spath, "STRU"))
    for j in filelist:
        file = Path(ppath, j)
        shutil.copy(file, spath)