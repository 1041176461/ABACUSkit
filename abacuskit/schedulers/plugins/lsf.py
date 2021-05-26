'''
Date: 2021-02-20 22:06:03
LastEditors: jiyuyang
LastEditTime: 2021-04-23 16:47:25
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from abacuskit.schedulers.datastructure import JobResource
from abacuskit.schedulers.schedulers import Scheduler

class LsfJobResource(JobResource):
    """
    An implementation of JobResource for LSF, that supports
    the OPTIONAL specification of a parallel environment (a string) + the total
    number of processors.

    'parallel_env' should contain a string of the form
    "host1 host2! hostgroupA! host3 host4" where the "!" symbol indicates the
    first execution host candidates. Other hosts are added only if the number of
    processors asked is more than those of the first execution host.
    See https://www-01.ibm.com/support/knowledgecenter/SSETD4_9.1.2/lsf_command_ref/bsub.1.dita?lang=en
    for more details about the parallel environment definition (the -m option of bsub).
    """

    _default_fields = (
        'parallel_env',
        'tot_num_mpiprocs',
    )

    @classmethod
    def validate_resources(cls, **kwargs):
        """Validate the resources against the job resource class of this scheduler."""
        resources = super().validate_resources(**kwargs)

        try:
            resources["parallel_env"] = str(kwargs.pop('parallel_env', ''))
        except (TypeError, ValueError):
            raise TypeError("When specified, 'parallel_env' must be a string")

        try:
            resources["tot_num_mpiprocs"] = int(kwargs.pop('tot_num_mpiprocs'))
        except (KeyError, ValueError):
            raise TypeError('tot_num_mpiprocs must be specified and must be an integer')

        if resources["tot_num_mpiprocs"] <= 0:
            raise ValueError('tot_num_mpiprocs must be >= 1')

        return resources

    def __init__(self, **kwargs):
        """Initialize the job resources from the passed arguments.

        :raises ValueError: if the resources are invalid or incomplete
        """
        self.resources = self.validate_resources(**kwargs)

class LsfScheduler(Scheduler):
    """
    Support for the IBM LSF scheduler
    'https://www-01.ibm.com/support/knowledgecenter/SSETD4_9.1.2/lsf_welcome.html'
    """
    def _get_submit_script_header(self, job_tmpl):
        """
        Return the submit script header. See the following manual
        https://www-01.ibm.com/support/knowledgecenter/SSETD4_9.1.2/lsf_command_ref/bsub.1.dita?lang=en
        for more details about the possible options to bsub, in particular for
        the parallel environment definition (with the -m option).
        """
        
        import string
        import re

        empty_line = ''

        lines = []
        if job_tmpl.submit_as_hold:
            lines.append('#BSUB -H')

        if job_tmpl.rerunnable:
            lines.append('#BSUB -r')
        else:
            lines.append('#BSUB -rn')

        if job_tmpl.email:
            lines.append(f'#BSUB -u {job_tmpl.email}')

        if job_tmpl.email_on_started:
            lines.append('#BSUB -B')
        if job_tmpl.email_on_terminated:
            lines.append('#BSUB -N')

        if job_tmpl.job_name:
            job_title = re.sub(r'[^a-zA-Z0-9_.-]+', '', job_tmpl.job_name)
            if not job_title or (job_title[0] not in string.ascii_letters + string.digits):
                job_title = f'j{job_title}'
            job_title = job_title[:128]
            lines.append(f'#BSUB -J "{job_title}"')

        if not job_tmpl.import_sys_environment:
            print('LSF scheduler cannot ignore the user environment', flush=True)

        if job_tmpl.sched_output_path:
            lines.append(f'#BSUB -o {job_tmpl.sched_output_path}')

        sched_error_path = getattr(job_tmpl, 'sched_error_path', None)
        if job_tmpl.sched_join_files:
            sched_error_path = f'{job_tmpl.sched_output_path}_'
            print(
                'LSF scheduler does not support joining '
                'the standard output and standard error '
                'files; std error file assigned instead '
                'to the file {}'.format(sched_error_path), flush=True
            )

        if sched_error_path:
            lines.append(f'#BSUB -e {job_tmpl.sched_error_path}')

        if job_tmpl.queue_name:
            lines.append(f'#BSUB -q {job_tmpl.queue_name}')

        if job_tmpl.priority:
            # Specifies user-assigned job priority that orders all jobs
            # (from all users) in a queue. Valid values for priority
            # are any integers between 1 and MAX_USER_PRIORITY
            # (configured in lsb.params, displayed by "bparams -l").
            # Jobs are scheduled based first on their queue priority first, then
            # job priority, and lastly in first-come first-served order.
            lines.append(f'#BSUB -sp {job_tmpl.priority}')

        if not job_tmpl.job_resource:
            raise ValueError('Job resources (as the tot_num_mpiprocs) are required for the LSF scheduler plugin')

        lines.append(f'#BSUB -n {job_tmpl.job_resource["tot_num_mpiprocs"]}')
        # Note:  make sure that PARALLEL_SCHED_BY_SLOT=Y is NOT
        # defined in lsb.params (you can check with the output of bparams -l).
        # Note: the -n option of bsub can also contain a maximum number of
        # procs to be used
        if job_tmpl.job_resource["parallel_env"]:
            lines.append(f'#BSUB -m "{job_tmpl.job_resource["parallel_env"]}"')

        if job_tmpl.max_wallclock_seconds is not None:
            # ABS_RUNLIMIT=Y should be set, in lsb.params (check with bparams -l)
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
            # The double negation results in the ceiling rather than the floor
            # of the division
            minutes = -(-(tot_secs % 3600) // 60)
            lines.append(f'#BSUB -W {hours:02d}:{minutes:02d}')

        if job_tmpl.max_memory_kb:
            try:
                virtual_memory_kb = int(job_tmpl.max_memory_kb)
                if virtual_memory_kb <= 0:
                    raise ValueError
            except ValueError:
                raise ValueError(
                    'max_memory_kb must be '
                    "a positive integer (in kB)! It is instead '{}'"
                    ''.format((job_tmpl.max_memory_kb))
                )
            # The -M option sets a per-process (soft) memory limit for all the
            # processes that belong to this job
            lines.append(f'#BSUB -M {virtual_memory_kb}')

        if job_tmpl.custom_scheduler_commands:
            lines.append(job_tmpl.custom_scheduler_commands)

        if job_tmpl.job_environment:
            lines.append(empty_line)
            lines.append('# ENVIRONMENT VARIABLES BEGIN ###')
            if not isinstance(job_tmpl.job_environment, dict):
                raise ValueError('If you provide job_environment, it must be a dictionary')
            for key, value in job_tmpl.job_environment.items():
                lines.append(f'export {key.strip()}={value}')
            lines.append('# ENVIRONMENT VARIABLES END  ###')
            lines.append(empty_line)

        lines.append(empty_line)

        # The following seems to be the only way to copy the input files
        # to the node where the computation are actually launched (the
        # -f option of bsub that does not always work...)
        # (need to add the line "#BSUB -outdir PATH_TO_REMOTE_DIRECTORY")
        # IMPORTANT! the -z is needed, because if LSB_OUTDIR is not defined,
        # you would do 'cp -R /* .' basically copying ALL FILES in your
        # computer (including mounted partitions) in the current dir!!
        lines.append("""
if [ ! -z "$LSB_OUTDIR" ]
then
  cp -R "$LSB_OUTDIR"/* .
fi
""")

        return '\n'.join(lines)

    def _get_submit_command(self, submit_script):
        """
        Return the string to execute to submit a given script.
        """
        submit_command = f'bsub < {submit_script}'

        return submit_command
