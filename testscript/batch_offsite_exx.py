from pyautotest.schedulers.data import Code
from pyautotest.calculations.plugins.exx import SetDimers, OptABFs, EXX
from pyautotest.utils.IO import write_input, read_kpt, read_stru

ntype = 4
ecut = 100
stru = read_stru(ntype, "STRU")
Nu = [6, 5, 4, 3, 2, 1]
kpt = read_kpt("KPT")
input_dict = {
                "ntype"         :   ntype,
                "ecutwfc"       :   ecut,
                "nbands"        :   200,
                "dr2"           :   1e-6,
                "niter"         :   400,
                "basis_type"    :   "lcao",
                "ks_solver"     :   "genelpa",
                "smearing"      :   "gaussian",
                "sigma"         :   0.0005,
                "mixing_type"   :   "pulay",
                "mixing_beta"   :   0.4,
                "pseudo_dir"    :   "/home/jiyy/SG15_upf",
                "newdm"         :   1,
                "exx_hybrid_type":         "hse",
                "exx_opt_orb_ecut":        ecut,
                "exx_opt_orb_tolerence":   1e-12,
                "exx_ccp_rmesh_times":     1.5,
                "exx_dm_threshold":        1e-5,
                "exx_cauchy_threshold":    1e-9,
                "exx_schwarz_threshold":   1e-5,
                "exx_c_threshold":         1e-5,
                "exx_v_threshold":         1e-3,
                "exx_pca_threshold":       1e-4
}
write_input(input_dict)
code_name="/home/jiyy/ABACUS/ABACUS_git_exx/git210119/bin/ABACUS.mpi.1.0.0_exx_fixdm3bug_avx2"
dimer_num = 5
# dimer calculation
obj_d = SetDimers(input_dict, stru, kpt, Nu, dimer_num)
command = Code(code_name=code_name,
	cmdline_params=["-n 1"],
	stdout_name="job.log",
	stderr_name="job.err",
	withmpi="mpirun"		
)
scheduler = "torque"
obj_d.batch_run(command, scheduler, num_machines=1, num_mpiprocs_per_machine=8, queue_name="gold5120", max_wallclock_seconds=360000, import_sys_environment=True)

# dimer check
obj_d.double_run(command, scheduler, num_machines=1, num_mpiprocs_per_machine=12, queue_name="gold5120", max_wallclock_seconds=360000, import_sys_environment=True)

# calculations below need submit script
# optmize orbitals
obj_a = OptABFs(ecut, stru, Nu)
external_command = "python -u /home/jiyy/ABACUS/ABACUS_git_exx/tools/opt_orb_pytorch/main.py > job.log 2>job.err"
obj_a.calculate(external_command)

# exx calculation
obj_e = EXX(input_dict, stru, kpt, Nu, dimer_num)
command = Code(code_name=code_name,
	cmdline_params=["-n 1", "-env OMP_NUM_THREADS=28"],
	stdout_name="job.log",
	stderr_name="job.err",
	withmpi="mpirun"		
)
obj_e.calculate(command)