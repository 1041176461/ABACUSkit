##################################################################################################################################
# Case 1 (Only PBE)

# 1. PBE cell-relax for ground structure to avoid the presence of imaginary frequencies in step 5
# 2. manually set `ocp` to 1 and `ocp_set` in INPUT for excited state calculation
# 3. PBE relax with INPUT given in step 2 based on the optimized structure given in step 1
# 4. using command `phonopy -d --dim="1 1 1"` to generate a series of structures with different displacements based on
# optimized structure given from step 1 and calculate force of them using PBE
# 5. using Phonopy to calculate the force matrix and phonon spectrum
# 6. using PyPhotonics to calculate PL shape based on the structures from step 1 and 3, and `band.yaml` file from step 5

##################################################################################################################################

##################################################################################################################################
# Case 2 (PBE and HSE)

# 1. PBE cell-relax for ground structure to avoid the presence of imaginary frequencies in step 5
# 2. HSE scf for optimized ground structure given from step 1
# 3. manually set `ocp` to 1 and `ocp_set` in INPUT for excited state calculation
# 4. PBE relax with INPUT given in step 2 based on the optimized structure given in step 1
# 5. HSE scf for optimized excited structure given from step 4
# 6. using command `phonopy -d --dim="1 1 1"` to generate a series of structures with different displacements based on
# optimized structure given from step 1 and calculate force of them using HSE
# 7. using Phonopy to calculate the force matrix and phonon spectrum
# 8. using PyPhotonics to calculate PL shape based on the structures from step 1 and 4, and `band.yaml` file from step 7

##################################################################################################################################

##################################################################################################################################
# Notes

# 1. `pulay-kerker` with a small `mixing_beta` for charge mixing may help for fast convergence of excited state calculation
# 2. if excited state calculation cannot converge within `relax_nmax`
# 3. excite one electron: ocp_set = [O_UP-1]*1.0 1*0.0 1*1.0 [U_UP-1]*0.0 <repeated N_K times> [O_DN]*1.0 [U_DN]*0.0 <repeated N_K times>

##################################################################################################################################

from abacuskit.utils.IO import read_stru
from abacuskit.utils.tools import list_elem2str
from abacuskit.schedulers.data import Code
from abacuskit.calculations.plugins.scf import RELAX, CELL_RELAX, SCF
from abacuskit.calculations.structure import Kpt
from abacuskit.utils.script import set_scheduler
from abacuskit.core.convert import Convert

import re
import os
import shutil
import numpy as np
from glob import glob
from copy import deepcopy
from textwrap import dedent
from phonopy.structure.atoms import atom_data, symbol_map


class Pre:

    def __init__(self, input_dict) -> None:
        self.input_dict = input_dict

    def _set_cell_relax_input(self, force_thr_ev=0.01, stress_thr=1):
        new_input_dict = deepcopy(self.input_dict)
        new_input_dict["relax_nmax"] = 50
        new_input_dict["cal_force"] = 1
        new_input_dict["force_thr_ev"] = force_thr_ev  # in eV/A
        new_input_dict["cal_stress"] = 1
        new_input_dict["stress_thr"] = stress_thr  # in KBar
        new_input_dict["relax_method"] = "cg"
        new_input_dict["out_force"] = 1
        new_input_dict["out_stru"] = 1
        return new_input_dict

    def _set_relax_input(self, force_thr_ev=0.01):
        new_input_dict = deepcopy(self.input_dict)
        new_input_dict["relax_nmax"] = 50
        new_input_dict["cal_force"] = 1
        new_input_dict["force_thr_ev"] = force_thr_ev  # in eV/A
        new_input_dict["relax_method"] = "cg"
        new_input_dict["out_force"] = 1
        new_input_dict["out_stru"] = 1
        return new_input_dict

    def ground_scf(self, stru, kpt):
        new_input_dict = deepcopy(self.input_dict)
        new_input_dict['cal_force'] = 1
        new_input_dict['out_force'] = 1
        obj = SCF(new_input_dict, stru, kpt)
        obj._prepare()

    def ground_relax(self, stru, kpt):
        new_input_dict = self._set_relax_input()
        obj = RELAX(new_input_dict, stru, kpt)
        obj._prepare()

    def ground_cell_relax(self, stru, kpt):
        new_input_dict = self._set_cell_relax_input()
        obj = CELL_RELAX(new_input_dict, stru, kpt)
        obj._prepare()

    def excite_relax(self, ocp_set, stru, kpt):
        new_input_dict = self._set_relax_input()
        new_input_dict['ocp'] = 1
        new_input_dict['ocp_set'] = ocp_set
        obj = RELAX(new_input_dict, stru, kpt)
        obj._prepare()


class Cal(Pre):

    def __init__(self, input_dict, code_name) -> None:
        super().__init__(input_dict)
        self.code_name = code_name

    def cal_ground_relax(self, knum, mpi=28, openmp=1, num_machines=1, num_mpiprocs_per_machine=28, max_wallclock_seconds=86400, queue_name='gold5120'):
        stru = read_stru(self.input_dict["ntype"], "STRU")
        kpt = Kpt('Gamma', knum)
        os.mkdir('ground')
        os.chdir('ground')
        command = Code(code_name=self.code_name, cmdline_params=[
                       f"-n {mpi}", f"-env OMP_NUM_THREADS={openmp}"], stdout_name="job.log", stderr_name="job.err", withmpi="mpirun")
        submit_command = set_scheduler('torque', [command], num_machines=num_machines, num_mpiprocs_per_machine=num_mpiprocs_per_machine,
                                       run_mode='p', max_wallclock_seconds=max_wallclock_seconds, queue_name=queue_name)
        self.ground_relax(stru, kpt)
        os.system(submit_command)
        os.chdir('../')

    def cal_ground_scf(self, knum, mpi=28, openmp=1, num_machines=1, num_mpiprocs_per_machine=28, max_wallclock_seconds=86400, queue_name='gold5120'):
        stru = read_stru(self.input_dict["ntype"], "STRU")
        kpt = Kpt('Gamma', knum)
        os.mkdir('ground')
        os.chdir('ground')
        command = Code(code_name=self.code_name, cmdline_params=[
                       f"-n {mpi}", f"-env OMP_NUM_THREADS={openmp}"], stdout_name="job.log", stderr_name="job.err", withmpi="mpirun")
        submit_command = set_scheduler('torque', [command], num_machines=num_machines, num_mpiprocs_per_machine=num_mpiprocs_per_machine,
                                       run_mode='p', max_wallclock_seconds=max_wallclock_seconds, queue_name=queue_name)
        self.ground_scf(stru, kpt)
        os.system(submit_command)
        os.chdir('../')

    def cal_ground_cell_relax(self, knum, mpi=28, openmp=1, num_machines=1, num_mpiprocs_per_machine=28, max_wallclock_seconds=86400, queue_name='gold5120'):
        stru = read_stru(self.input_dict["ntype"], "STRU")
        kpt = Kpt('Gamma', knum)
        os.mkdir('ground')
        os.chdir('ground')
        command = Code(code_name=self.code_name, cmdline_params=[
                       f"-n {mpi}", f"-env OMP_NUM_THREADS={openmp}"], stdout_name="job.log", stderr_name="job.err", withmpi="mpirun")
        submit_command = set_scheduler('torque', [command], num_machines=num_machines, num_mpiprocs_per_machine=num_mpiprocs_per_machine,
                                       run_mode='p', max_wallclock_seconds=max_wallclock_seconds, queue_name=queue_name)
        self.ground_cell_relax(stru, kpt)
        os.system(submit_command)
        os.chdir('../')

    def param_for_ocp_set(self):
        """Only spin 1"""
        filename = glob('./ground/OUT*/running_*.log')
        if filename:
            with open(filename[0], 'r') as f:
                o_up = int(search_line(f, "occupied bands"))
                nbands = int(search_line(f, 'NBANDS'))
                u_up = nbands - o_up
                nk = int(search_line(f, 'nkstot'))
        else:
            raise FileNotFoundError(
                "Not found log file to extract parameters for ocp_set.")
        return nbands, o_up, u_up, nk

    def exite_string_ocp_set(self, o, u, nk, enum=1):
        return f"{o-enum}*1.0 {enum}*0.0 {enum}*1.0 {u-enum}*0.0 "*nk

    def ground_string_ocp_set(self, o, u, nk):
        return f"{o}*1.0 {u}*0.0 "*nk

    def cal_excite_relax(self, knum, ocp_set, mpi=28, openmp=1, num_machines=1, num_mpiprocs_per_machine=28, max_wallclock_seconds=86400, queue_name='gold5120'):
        strufile = extract_largest_number_file()
        stru = read_stru(self.input_dict["ntype"], glob(strufile))
        for elem in stru.orbitals.keys():
            stru.orbitals[elem] = os.path.basename(stru.orbitals[elem])
        kpt = Kpt('Gamma', knum)
        os.mkdir('excite')
        os.chdir('excite')
        command = Code(code_name=self.code_name, cmdline_params=[
                       f"-n {mpi}", f"-env OMP_NUM_THREADS={openmp}"], stdout_name="job.log", stderr_name="job.err", withmpi="mpirun")
        submit_command = set_scheduler('torque', [command], num_machines=num_machines, num_mpiprocs_per_machine=num_mpiprocs_per_machine,
                                       run_mode='p', max_wallclock_seconds=max_wallclock_seconds, queue_name=queue_name)
        self.excite_relax(ocp_set, stru, kpt)
        os.system(submit_command)
        os.chdir('../')

    def check_force(self, filename):
        count = 1
        with open(filename, 'r') as file:
            for line in file:
                if re.search(r"TOTAL ATOM NUMBER = [0-9]+", line):
                    natom = int(re.search("[0-9]+", line).group())
                    force = np.zeros((natom, 3))
                if re.search("TOTAL-FORCE \(eV/Angstrom\)", line):
                    for i in range(4):
                        file.readline()
                    for i in range(natom):
                        _, fx, fy, fz = file.readline().split()
                        force[i] = (float(fx), float(fy), float(fz))
                    print(f"Variance of force in step {count}: {np.var(force)}")
                    count += 1
        return np.var(force)

    def zpl(self, glog, elog):
        with open(glog, 'r') as f1, open(elog, 'r') as f2:
            gE = search_line(f1, 'E_KohnSham')
            eE = search_line(f2, 'E_KohnSham')
        return float(eE)-float(gE)

    def cal_phono(self, dim=[1, 1, 1], mpi=28, openmp=1, num_machines=1, num_mpiprocs_per_machine=28, max_wallclock_seconds=86400, queue_name='gold5120'):
        strufile = extract_largest_number_file()
        stru = read_stru(self.input_dict["ntype"], strufile)
        for elem in stru.orbitals.keys():
            stru.orbitals[elem] = os.path.basename(stru.orbitals[elem])
        self.input_dict['gamma_only'] = 1
        kpt = Kpt('Gamma', [1, 1, 1])
        command = Code(code_name=self.code_name, cmdline_params=[
                       f"-n {mpi}", f"-env OMP_NUM_THREADS={openmp}"], stdout_name="job.log", stderr_name="job.err", withmpi="mpirun")
        os.mkdir('phonon_ground')
        os.chdir('phonon_ground')
        stru.write_stru()
        with open("setting.conf", 'w') as f:
            f.write(dedent(f"""\
                 DIM = {' '.join(list_elem2str(dim))}
                 ATOM_NAME = {' '.join(stru.elements)}
                """
                           ))
        os.system('phonopy setting.conf --abacus -d')
        disp_list = glob('STRU-*')
        for i in disp_list:
            dst = 'disp-'+i.split('-')[-1]
            os.mkdir(dst)
            shutil.move(i, dst)
            os.chdir(dst)
            disp_stru = read_stru(self.input_dict["ntype"], i)
            submit_command = set_scheduler('torque', [command], num_machines=num_machines, num_mpiprocs_per_machine=num_mpiprocs_per_machine,
                                           run_mode='p', max_wallclock_seconds=max_wallclock_seconds, queue_name=queue_name)
            self.ground_scf(disp_stru, kpt)
            os.system(submit_command)
            os.chdir('../')
        os.chdir('../')

    def create_force_sets(self):
        os.chdir('phonon_ground')
        #logfiles = glob('disp-*/OUT*/running*.log')
        # logfiles.reverse()
        #logfiles = ' '.join(logfiles)
        #os.system(f'phonopy -f {logfiles}')
        os.system('phonopy -f ./disp-*/OUT*/running*.log')
        os.chdir('../')

    def cal_phono_band(self, path, labels, points=51):
        os.chdir('phonon_ground')
        stru = read_stru(self.input_dict["ntype"], 'STRU')
        with open("band.conf", 'w') as f:
            f.write(dedent(f"""\
                 ATOM_NAME = {' '.join(stru.elements)}
                 BAND = {convert_path_to_string(path)}
                 BAND_POINTS = {points}
                 BAND_LABELS = {' '.join(labels)}
                 BAND_CONNECTION = .TRUE.
                 FORCE_CONSTANTS= WRITE 
                """
                           ))
        os.system("phonopy band.conf -p -s")
        os.chdir('../')

    def cal_emission(self, num_mode):
        stru = read_stru(self.input_dict["ntype"], 'phonon_ground/STRU')
        mass = []
        for elem in stru.elements:
            mass += [atom_data[symbol_map[elem]][3]]*stru.numbers[elem]
        mass = np.array(mass)*1.660539040e-27 # in unit AMU
        os.mkdir('photon')
        shutil.copyfile(extract_largest_number_file(dst='./ground/OUT*'), 'photon/STRU_GS')
        shutil.copyfile(extract_largest_number_file(dst='./excite/OUT*'), 'photon/STRU_ES')
        shutil.copy2('phonon_ground/band.yaml', 'photon/')
        os.chdir('photon')
        Convert.STRU_to_POSCAR('STRU_GS', 'CONTCAR_GS', self.input_dict["ntype"])
        Convert.STRU_to_POSCAR('STRU_ES', 'CONTCAR_ES', self.input_dict["ntype"])

        from pyphotonics.photoluminescence import Photoluminescence
        p = Photoluminescence('./', 'CONTCAR_GS', 'CONTCAR_ES', num_mode, "phonopy", mass, 1000, shift_vector=[0.0, 0.0, 0.0])

        import matplotlib.pyplot as plt
        print("Delta_R=", p.Delta_R)
        print("Delta_Q=", p.Delta_Q)
        print("HuangRhyes=", p.HuangRhyes)

        plt.figure(figsize=(10, 10))
        plt.plot(p.S_omega)
        plt.ylabel('$S(\hbar\omega)$')
        plt.xlabel('Phonon energy (meV)')
        plt.xlim(0, 200)
        # plt.ylim(0, 0.01)
        plt.savefig('S_omega', bbox_inches='tight')
        p.write_S('S')

        A, I = p.PL(2, 2, 1.95)
        plt.figure(figsize=(10, 10))
        plt.plot(I.__abs__())
        plt.ylabel('$I(\hbar\omega)$')
        plt.xlabel('Photon energy (eV)')
        plt.xlim(1200, 2000)
        x_values, labels = plt.xticks()
        labels = [float(x)/p.resolution for x in x_values]
        plt.xticks(x_values, labels)
        plt.ylim(0, 600)
        plt.savefig('I', bbox_inches='tight')

        os.chdir('../')


def convert_path_to_string(path=[[]]):
    line = []
    for i in path:
        line.append(' '.join(list_elem2str(i)))
    return '\t'.join(line)


def extract_ion_number_from_file(s):
    return int(re.search("(ION)([0-9]+)", s).group(2))


def extract_largest_number_file(dst='./ground/OUT*'):
    filelist = glob(os.path.join(dst, 'STRU_ION*_D'))
    dstname = os.path.dirname(filelist[-1])
    numlist = list(map(extract_ion_number_from_file, filelist))
    numlist.sort()
    return os.path.join(dstname, f'STRU_ION{numlist[-1]}_D')


def search_line(f, name):
    for line in f:
        if re.search(name, line):
            res = line.split()[-1]
    f.seek(0)
    return res


if __name__ == "__main__":
    code_name = "/home/jiyy/ABACUS/abacus-develop-2.2.0/bin/ABACUS.mpi"
    input_dict = {
        "ntype":   2,
        "ecutwfc":   100,
        "nbands":   220,
        "scf_nmax":   400,
        "scf_thr":   1e-6,
        "basis_type":   "lcao",
        "ks_solver":   "genelpa",
        "smearing_method":   "gaussian",
        "smearing_sigma":   0.0005,
        "mixing_type":   "pulay-kerker",
        "mixing_beta":   0.4,
        "pseudo_dir":   "/home/jiyy/potential/SG15/SG15_upf-soc",
        "orbital_dir":   "/home/jiyy/GD_orb/dpsi-T-S",
        "out_chg":   1,
        "nspin":   2,
        "nelec":   254,
        "gamma_only":   1
    }
    num_machines = 1
    num_mpiprocs_per_machine = 8
    max_wallclock_seconds = 86400
    mpi = 1
    openmp = 8
    knum = [1, 1, 1]
    dim = [1, 1, 1]

    cal = Cal(input_dict, code_name)
    # cal.cal_ground_relax(
    #    knum, mpi, openmp, num_machines, num_mpiprocs_per_machine, max_wallclock_seconds)
    #nbands, o_up, u_up, nk = cal.param_for_ocp_set()
    #print(cal.exite_string_ocp_set(o_up, u_up, nk))
    # cal.cal_excite_relax(
    #    knum, ocp_set, mpi, openmp, num_machines, num_mpiprocs_per_machine, max_wallclock_seconds)
    #cal.cal_phono(dim, mpi, openmp, num_machines,
    #              num_mpiprocs_per_machine, max_wallclock_seconds)
    cal.cal_emission(num_mode=189)
