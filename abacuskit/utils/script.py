'''
Date: 2021-03-16 20:34:00
LastEditors: jiyuyang
LastEditTime: 2021-04-28 18:06:07
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from abacuskit.schedulers.data import Code, CodeRunMode, JobDefaultFields

import os
import json

def set_scheduler(
				scheduler,
				codes_info,
				num_machines=None,
				num_mpiprocs_per_machine=None,
				num_cores_per_machine=None,
				num_cores_per_mpiproc=None,
				parallel_env='',
				tot_num_mpiprocs=None,
				shebang=None,
                submit_as_hold=None,
                rerunnable=None,
                job_environment=None,
                working_directory=None,
                email=None,
                email_on_started=None,
                email_on_terminated=None,
                job_name=None,
                sched_output_path='job.log',
                sched_error_path='job.err',
                sched_join_files=False,
                queue_name=None,
                account=None,
                qos=None,
                priority=None,
                max_memory_kb=None,
                max_wallclock_seconds=86400,
                custom_scheduler_commands=None,
                prepend_text=None,
                append_text=None,
                import_sys_environment=True,
                run_mode='s',
				submit_script='sub.sh'
				):
	""" 
	Write the subnit script as a file.
	
	scheduler: choose a scheduler ('pbspro')
	codes_info: a list of which the (multiple) codes have to be executed.

	Specify the job resources:
	For PBSpro:
		num_machines: how many machines (or nodes) should be used
        num_mpiprocs_per_machine: how many MPI procs should be used on each machine (or node).
        num_cores_per_machine: how many cores should be used on each machine (or node).
		num_cores_per_mpiproc:

	For Torque:
		num_machines: how many machines (or nodes) should be used
        num_cores_per_machine/num_mpiprocs_per_machine: how many MPI procs should be used on each machine (or node).
		num_cores_per_mpiproc:

	For SLURM:
		num_machines: how many machines (or nodes) should be used
        num_mpiprocs_per_machine: how many MPI procs should be used on each machine (or node).	
        num_cores_per_mpiproc: (optional) how many cores should be used on MPI procs.
		num_cores_per_mpiproc:

	For LSF:
		tot_num_mpiprocs:
		parallel_env:

	For SGE:
		tot_num_mpiprocs:
		parallel_env:

	Others:
		shebang: the first line of the submission script. Default: None
        submit_as_hold: if set, the job will be in a 'hold' status right after the submission. Default: None
        rerunnable: if the job is rerunnable (boolean). Default: None
        job_environment: a dictionary with environment variables to set before the execution of the code. Default: {'OMP_NUM_THREADS':1}
        working_directory: the working directory for this job. Default: None
        email: an email address for sending emails on job events. Default: None
        email_on_started: if True, ask the scheduler to send an email when the job starts. Default: None
        email_on_terminated: if True, ask the scheduler to send an email when the job ends. Default: None
        job_name: the name of this job. Default: None
        sched_output_path: a (relative) file name for the stdout of this job. Default: 'job.log'
        sched_error_path: a (relative) file name for the stdout of this job. Default: 'job.err'
        sched_join_files: if True, write both stdout and stderr on the same file (the one specified for stdout). Default: False
        queue_name: the name of the scheduler queue (sometimes also called partition), on which the job will be submitted. Default: None
        account: the name of the scheduler account (sometimes also called projectid), on which the job will be submitted. Default: None
        qos: the quality of service of the scheduler account, on which the job will be submitted. Default: None
        priority: a priority for this job. Default: None
        max_memory_kb: the maximum amount of memory the job is allowed to allocate ON EACH NODE, in kilobytes. Default: None
        max_wallclock_seconds: the maximum wall clock time that all processes of a job are allowed to exist, in seconds. Default: 86400
        custom_scheduler_commands: a string that will be inserted right after the last scheduler command, and before any other non-scheduler command. Default: "ulimit -s unlimited"
        prepend_text: a (possibly multi-line) string to be inserted in the scheduler script before the main execution line. Default: None
        append_text: a (possibly multi-line) string to be inserted in the scheduler script after the main execution line. Default: None
        import_sys_environment: import the system environment variables. Default: True
        run_mode: it contains the information necessary to run a single code. 's': serial, 'p': parallel. Default: 'p'
	"""

	if run_mode == 'p':
		codes_run_mode = CodeRunMode.PARALLEL
	elif run_mode == 's':
		codes_run_mode = CodeRunMode.SERIAL
		
	if scheduler == 'pbspro':
		from abacuskit.schedulers.plugins.pbsbaseclasses import PbsJobResource
		from abacuskit.schedulers.plugins.pbspro import PbsproScheduler
		job_resource = PbsJobResource(num_machines=num_machines, num_mpiprocs_per_machine=num_mpiprocs_per_machine, num_cores_per_machine=num_cores_per_machine, num_cores_per_mpiproc=num_cores_per_mpiproc).resources
		sche = PbsproScheduler()

	elif scheduler == 'torque':
		from abacuskit.schedulers.plugins.pbsbaseclasses import PbsJobResource
		from abacuskit.schedulers.plugins.torque import TorqueScheduler
		job_resource = PbsJobResource(num_machines=num_machines, num_mpiprocs_per_machine=num_mpiprocs_per_machine, num_cores_per_machine=num_cores_per_machine, num_cores_per_mpiproc=num_cores_per_mpiproc).resources
		sche = TorqueScheduler()

	elif scheduler == 'slurm':
		from abacuskit.schedulers.plugins.slurm import SlurmJobResource, SlurmScheduler
		job_resource = SlurmJobResource(num_machines=num_machines, num_mpiprocs_per_machine=num_mpiprocs_per_machine, num_cores_per_machine=num_cores_per_machine, num_cores_per_mpiproc=num_cores_per_mpiproc).resources
		sche = SlurmScheduler()

	elif scheduler == 'lsf':
		from abacuskit.schedulers.plugins.lsf import LsfJobResource, LsfScheduler
		job_resource = LsfJobResource(parallel_env=parallel_env, tot_num_mpiprocs=tot_num_mpiprocs).resources
		sche = LsfScheduler()

	elif scheduler == 'sge':
		from abacuskit.schedulers.plugins.sge import SgeJobResource, SgeScheduler
		job_resource = SgeJobResource(parallel_env=parallel_env, tot_num_mpiprocs=tot_num_mpiprocs).resources
		sche = SgeScheduler()

	elif scheduler == 'direct':
		from abacuskit.schedulers.plugins.direct import DirectJobResource, DirectScheduler
		job_resource = DirectJobResource(num_machines=num_machines, num_mpiprocs_per_machine=num_mpiprocs_per_machine, tot_num_mpiprocs=tot_num_mpiprocs).resources
		sche = DirectScheduler()

	else:
		raise KeyError(f"It does not support scheduler `{scheduler}`")	

	job_tml = JobDefaultFields(
		shebang=shebang,
		submit_as_hold=submit_as_hold,
		rerunnable=rerunnable,
		job_environment=job_environment,
		working_directory=working_directory,
		email=email,
		email_on_started=email_on_started,
		email_on_terminated=email_on_terminated,
		job_name=job_name,
		sched_output_path=sched_output_path,
		sched_error_path=sched_error_path,
		sched_join_files=sched_join_files,
		queue_name=queue_name,
		account=account,
		qos=qos,
		job_resource=job_resource,
		priority=priority,
		max_memory_kb=max_memory_kb,
		max_wallclock_seconds=max_wallclock_seconds,
		custom_scheduler_commands=custom_scheduler_commands,
		prepend_text=prepend_text,
		append_text=append_text,
		import_sys_environment=import_sys_environment,
		codes_run_mode=codes_run_mode,
		codes_info=codes_info
	)
	
	out_string = sche.get_submit_script(job_tml)
	with open(submit_script, 'w') as file:
		file.write(out_string)

	return sche._get_submit_command(submit_script)

def submit_script(filename, line, **kwargs):
    """Run with submit script
    
    :params filename: string of input file name
    :params line: command line in submit script
    """
    
    with open(filename, 'r') as file:
        text = json.load(file)

    job_resource = text["job_resource"]
    script_params = text["script_params"]
    job_environment = text.pop("job_environment", None) #{"MKL_THREADING_LAYER":"GNU"} MKL_THREADING_LAYER=INTEL is incompatible with libgomp.so.1 library. See https://github.com/pytorch/pytorch/issues/37377
    cmdline_params = [line]
    codes_info = [Code(cmdline_params=cmdline_params,
						stdin_name=kwargs.pop("stdin_name", None),
						stdout_name=kwargs.pop("stdout_name", None),
						stderr_name=kwargs.pop("stderr_name", None),
						join_files=kwargs.pop("join_files", False),
						withmpi=kwargs.pop("withmpi", None),
						code_uuid=kwargs.pop("code_uuid", None))]
    submit_command = set_scheduler(codes_info=codes_info, job_environment=job_environment, **script_params, **job_resource)
    os.system(submit_command)
