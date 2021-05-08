'''
Date: 2021-03-18 20:12:03
LastEditors: jiyuyang
LastEditTime: 2021-04-28 19:00:25
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from enum import IntEnum

from numpy import savez_compressed

class Code:
    def __init__(self,
                code_name='',
                cmdline_params=[],
                stdin_name=None,
                stdout_name=None,
                stderr_name=None,
                join_files=None,
                withmpi='mpirun',
                code_uuid=None
    ):
        """
        params: code_name: string of code_name
        params: cmdline_params: a list of strings with the command line arguments of the program to run. For example: mpirun cmdline_params[0] cmdline_params[1] ... < stdin > stdout. 
        params: stdin_name: (optional) the name of the standard input file.  code.x < stdin_name
        params: stdout_name: (optional) the name of the standard output file. code.x ... > stdout_name
        params: stderr_name: (optional) a string, the name of the error file of the code.
        params: join_files: If join_files=True, code.x ... > stdout_name 2>&1, if join_files=False, code.x ... > stdout_name 2> stderr_name
        params: withmpi: if not None, executes the code with `withmpi` value (mpirun or another MPI installed on the remote computer)
        params: code_uuid: the uuid of the code
        """

        self.code_name = code_name
        self.cmdline_params = cmdline_params     
        self.stdin_name = stdin_name            
        self.stdout_name = stdout_name         
        self.stderr_name = stderr_name          
        self.join_files = join_files             
        self.withmpi = withmpi                   
        self.code_uuid = code_uuid

    def run_line(self):
        """Return run lines """

        command_to_exec = f'{self.withmpi} '+' '.join(self.cmdline_params)+f' {self.code_name}' if self.withmpi else ' '.join(self.cmdline_params)
        stdin_str = f'< {self.stdin_name}' if self.stdin_name else ''
        stdout_str = f'> {self.stdout_name}' if self.stdout_name else ''
        stderr_str = '2>&1' if self.join_files else f'2> {self.stderr_name}' if self.stderr_name else ''
        output_string = f'{command_to_exec} {stdin_str} {stdout_str} {stderr_str}'
        
        return output_string

class CodeRunMode(IntEnum):
    """Enum to indicate the way the codes of a calculation should be run.

    For PARALLEL, the codes for a given calculation will be run in parallel by running them in the background::
        
        code1.x &
        
        code2.x &
    
    For the SERIAL option, codes will be executed sequentially by running for example the following::
        
        code1.x
        code2.x
        
    """

    SERIAL = 0
    PARALLEL = 1

class JobDefaultFields:
    def __init__(self,
                shebang=None,
                submit_as_hold=None,
                rerunnable=None,
                job_environment=None,
                working_directory=None,
                email=None,
                email_on_started=None,
                email_on_terminated=None,
                job_name=None,
                sched_output_path=None,
                sched_error_path=None,
                sched_join_files=None,
                queue_name=None,
                account=None,
                qos=None,
                job_resource=None,
                priority=None,
                max_memory_kb=None,
                max_wallclock_seconds=None,
                custom_scheduler_commands=None,
                prepend_text=None,
                append_text=None,
                import_sys_environment=None,
                codes_run_mode=None,
                codes_info=None
                ):

        self.shebang = shebang                                      # The first line of the submission script.
        self.submit_as_hold = submit_as_hold                        # if set, the job will be in a 'hold' status right after the submission.
        self.rerunnable = rerunnable                                # if the job is rerunnable (boolean).
        self.job_environment = job_environment                      # a dictionary with environment variables to set before the execution of the code
        self.working_directory = working_directory                  # the working directory for this job.
        self.email = email                                          # an email address for sending emails on job events.
        self.email_on_started = email_on_started                    # if True, ask the scheduler to send an email when the job starts.
        self.email_on_terminated = email_on_terminated              # if True, ask the scheduler to send an email when the job ends.
        self.job_name = job_name                                    # the name of this job.
        self.sched_output_path = sched_output_path                  # a (relative) file name for the stdout of this job
        self.sched_error_path = sched_error_path                    # a (relative) file name for the stdout of this job
        self.sched_join_files = sched_join_files                    # if True, write both stdout and stderr on the same file (the one specified for stdout)
        self.queue_name = queue_name                                # the name of the scheduler queue (sometimes also called partition), on which the job will be submitted.
        self.account = account                                      # the name of the scheduler account (sometimes also called projectid), on which the job will be submitted.
        self.qos = qos                                              # the quality of service of the scheduler account, on which the job will be submitted.
        self.job_resource = job_resource                            # how many nodes and cpus it should use.
        self.priority = priority                                    # a priority for this job.
        self.max_memory_kb = max_memory_kb                          # The maximum amount of memory the job is allowed to allocate ON EACH NODE, in kilobytes.
        self.max_wallclock_seconds = max_wallclock_seconds          # The maximum wall clock time that all processes of a job are allowed to exist, in seconds.
        self.custom_scheduler_commands = custom_scheduler_commands  # a string that will be inserted right after the last scheduler command, and before any other non-scheduler command.
        self.prepend_text = prepend_text                            # a (possibly multi-line) string to be inserted in the scheduler script before the main execution line.
        self.append_text = append_text                              # a (possibly multi-line) string to be inserted in the scheduler script after the main execution line.
        self.import_sys_environment = import_sys_environment        # import the system environment variables.
        self.codes_run_mode = codes_run_mode                        # it contains the information necessary to run a single code.
        self.codes_info = codes_info                                # sets which the (multiple) codes have to be executed.