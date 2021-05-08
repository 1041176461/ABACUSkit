from pyautotest.calculations.structure import read_stru, read_kpt
from pyautotest.schedulers.data import Code
from pyautotest.calculations.plugins.exx import SetDimers
from pyautotest.calculations.baseclass import ABACUSCalculation
ntype = 4
stru = read_stru(ntype, "STRU")
kpt = read_kpt("KPT")
input_dict = {
                "ntype"         :   ntype,
                "ecutwfc"       :   100,
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
                "exx_opt_orb_ecut":        100,
                "exx_opt_orb_tolerence":   1e-12,
                "exx_ccp_rmesh_times":     1.5,
                "exx_dm_threshold":        1e-5,
                "exx_cauchy_threshold":    1e-9,
                "exx_schwarz_threshold":   1e-5,
                "exx_c_threshold":         1e-5,
                "exx_v_threshold":         1e-3,
                "exx_pca_threshold":       1e-4
}
with open("INPUT", 'w') as file:
    file.write(ABACUSCalculation().get_input_line(input_dict))
Nu = [6, 5, 4, 3, 2, 1]
dimer_num = 5
obj = SetDimers(input_dict, stru, kpt, Nu, dimer_num)
command = Code(code_name="/home/jiyy/ABACUS/ABACUS_git_exx/git210119/bin/ABACUS.mpi.1.0.0_exx_fixdm3bug_avx2",
	cmdline_params=["-n 1"],
	stdout_name="job.log",
	stderr_name="job.err",
	withmpi="mpirun"		
)
scheduler = "torque"
obj.batch_run(command, scheduler, num_machines=1, num_mpiprocs_per_machine=16, queue_name="gold5120", max_wallclock_seconds=360000, import_sys_environment=True)
# obj.double_run(command, scheduler, **kwargs)