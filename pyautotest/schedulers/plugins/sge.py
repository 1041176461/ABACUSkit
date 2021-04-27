'''
Date: 2021-02-21 13:03:33
LastEditors: jiyuyang
LastEditTime: 2021-04-23 16:47:37
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.schedulers.schedulers import Scheduler
from pyautotest.schedulers.datastructure import ParEnvJobResource

class SgeJobResource(ParEnvJobResource):
    pass

class SgeScheduler(Scheduler):
    """
    Support for the Sun Grid Engine scheduler and its variants/forks (Son of Grid Engine, Oracle Grid Engine, ...)
    """

    def _get_submit_script_header(self, job_tmpl):
        """
        Return the submit script header
        """
        import re
        import string

        empty_line = ''

        lines = []

        # SGE provides flags for wd and cwd
        if job_tmpl.working_directory:
            lines.append(f'#$ -wd {job_tmpl.working_directory}')
        else:
            lines.append('#$ -cwd')

        # Enforce bash shell
        lines.append('#$ -S /bin/bash')

        if job_tmpl.submit_as_hold:
            lines.append(f'#$ -h {job_tmpl.submit_as_hold}')

        if job_tmpl.rerunnable:
            lines.append(f'#$ -r {job_tmpl.rerunnable}')

        if job_tmpl.email:
            lines.append(f'#$ -M {job_tmpl.email}')

        email_events = ''
        if job_tmpl.email_on_started:
            email_events += 'b'
        if job_tmpl.email_on_terminated:
            email_events += 'ea'
        if email_events:
            lines.append(f'#$ -m {email_events}')
            if not job_tmpl.email:
                print(
                    'Email triggers provided to SGE script for job,'
                    'but no email field set; will send emails to '
                    'the job owner as set in the scheduler', flush=True
                )
        else:
            lines.append('#$ -m n')

        # From the qsub man page:
        # "The name may be any arbitrary alphanumeric ASCII string, but
        # may  not contain  "\n", "\t", "\r", "/", ":", "@", "\", "*",
        # or "?"."
        if job_tmpl.job_name:
            job_title = re.sub(r'[^a-zA-Z0-9_.-]+', '', job_tmpl.job_name)
            if not job_title or (job_title[0] not in string.ascii_letters):
                job_title = f'j{job_title}'

            lines.append(f'#$ -N {job_tmpl.job_name}')

        if job_tmpl.import_sys_environment:
            lines.append('#$ -V')

        if job_tmpl.sched_output_path:
            lines.append(f'#$ -o {job_tmpl.sched_output_path}')

        if job_tmpl.sched_join_files:
            # from qsub man page:
            # 'y': Standard error and standard output are merged  into
            #       standard output
            # 'n' : Standard error and standard output are not merged (default)
            lines.append('#$ -j y')
            if job_tmpl.sched_error_path:
                print(
                    'sched_join_files is True, but sched_error_path is set in '
                    'SGE script; ignoring sched_error_path', flush=True
                )
        else:
            if job_tmpl.sched_error_path:
                lines.append(f'#$ -e {job_tmpl.sched_error_path}')

        if job_tmpl.queue_name:
            lines.append(f'#$ -q {job_tmpl.queue_name}')

        if job_tmpl.account:
            lines.append(f'#$ -P {job_tmpl.account}')

        if job_tmpl.priority:
            # Priority of the job.  Format: host-dependent integer.  Default:
            # zero.   Range:  [-1023,  +1024].  Sets job's Priority
            # attribute to priority.
            lines.append(f'#$ -p {job_tmpl.priority}')

        if not job_tmpl.job_resource:
            raise ValueError('Job resources (as the tot_num_mpiprocs) are required for the SGE scheduler plugin')
        # Setting up the parallel environment
        lines.append(f'#$ -pe {str(job_tmpl.job_resource["parallel_env"])} {int(job_tmpl.job_resource["tot_num_mpiprocs"])}')

        if job_tmpl.max_wallclock_seconds is not None:
            try:
                tot_secs = int(job_tmpl.max_wallclock_seconds)
                if tot_secs <= 0:
                    raise ValueError
            except ValueError:
                raise ValueError(
                    'max_wallclock_seconds must be '
                    "a positive integer (in seconds)! It is instead '{}'"
                    ''.format((job_tmpl.max_wallclock_seconds))
                )
            hours = tot_secs // 3600
            tot_minutes = tot_secs % 3600
            minutes = tot_minutes // 60
            seconds = tot_minutes % 60
            lines.append(f'#$ -l h_rt={hours:02d}:{minutes:02d}:{seconds:02d}')

        if job_tmpl.custom_scheduler_commands:
            lines.append(job_tmpl.custom_scheduler_commands)

        if job_tmpl.job_environment:
            lines.append(empty_line)
            lines.append('# ENVIRONMENT VARIABLES BEGIN ###')
            if not isinstance(job_tmpl.job_environment, dict):
                raise ValueError('If you provide job_environment, it must be a dictionary')
            for key, value in job_tmpl.job_environment.items():
                lines.append(f'export {key.strip()}={value}')
            lines.append('# ENVIRONMENT VARIABLES  END  ###')
            lines.append(empty_line)

        return '\n'.join(lines)

    def _get_submit_command(self, submit_script):
        """
        Return the string to execute to submit a given script.
        """

        submit_command = f'qsub -terse {submit_script}'

        return submit_command
