'''
Date: 2021-02-12 11:57:16
LastEditors: jiyuyang
LastEditTime: 2021-04-23 16:46:29
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from abacuskit.schedulers.datastructure import NodeNumberJobResource
from abacuskit.schedulers.schedulers import Scheduler

class PbsJobResource(NodeNumberJobResource):
    """Class for PBS job resources."""

    @classmethod
    def validate_resources(cls, **kwargs):
        """Validate the resources against the job resource class of this scheduler.

        This extends the base class validator and calculates the `num_cores_per_machine` fields to pass to PBSlike
        schedulers. Checks that `num_cores_per_machine` is a multiple of `num_cores_per_mpiproc` and/or
        `num_mpiprocs_per_machine`.

        :param kwargs: dictionary of values to define the job resources
        :return: defaultdict with the parsed parameters populated
        :raises ValueError: if the resources are invalid or incomplete
        """
        resources = super().validate_resources(**kwargs)

        if resources['num_cores_per_machine'] is not None and resources['num_cores_per_mpiproc'] is not None:
            if resources['num_cores_per_machine'] != resources['num_cores_per_mpiproc'] * resources['num_mpiprocs_per_machine']:
                raise ValueError(
                    '`num_cores_per_machine` must be equal to `num_cores_per_mpiproc * num_mpiprocs_per_machine` and in'
                    ' particular it should be a multiple of `num_cores_per_mpiproc` and/or `num_mpiprocs_per_machine`'
                )

        elif resources['num_cores_per_mpiproc'] is not None:
            if resources['num_cores_per_mpiproc'] < 1:
                raise ValueError('num_cores_per_mpiproc must be greater than or equal to one.')
            
            # In this plugin we never used num_cores_per_mpiproc so if it is not defined it is OK.
            resources['num_cores_per_machine'] = (resources['num_cores_per_mpiproc'] * resources['num_mpiprocs_per_machine'])

        return resources

class PbsBaseClass(Scheduler):
    """Base class with support for the PBSPro scheduler
    (http://www.pbsworks.com/) and for PBS and Torque
    (http://www.adaptivecomputing.com/products/open-source/torque/).
    """

    def _get_resource_lines(
        self, num_machines, num_mpiprocs_per_machine, num_cores_per_machine, max_memory_kb, max_wallclock_seconds
    ):
        """
        Return a set a list of lines (possibly empty) with the header
        lines relative to:

        * num_machines
        * num_mpiprocs_per_machine
        * num_cores_per_machine
        * max_memory_kb
        * max_wallclock_seconds
        
        This is done in an external function because it may change in
        different subclasses.
        """
        raise NotImplementedError('Implement the _get_resource_lines in each subclass!')

    def _get_submit_script_header(self, job_tmpl):
        """
        Return the submit script header, using the parameters from the job_tmpl.
        """
        import re
        import string
        empty_line = ''

        lines = []
        if job_tmpl.submit_as_hold:
            lines.append('#PBE -h')
        
        if job_tmpl.rerunnable:
            lines.append('#PBE -r y')
        else:
            lines.append('#PBE -r n')
        
        if job_tmpl.email:
            lines.append(f'#PBE -M {job_tmpl.email}')

        email_events = ''
        if job_tmpl.email_on_started:
            email_events += 'b'
        if job_tmpl.email_on_terminated:
            email_events += 'ea'
        if email_events:
            lines.append(f'#PBS -m {email_events}')
            if not job_tmpl.email:
                print(
                    'Email triggers provided to PBSPro script for job,'
                    'but no email field set; will send emails to '
                    'the job owner as set in the scheduler', flush=True
                )
        else:
            lines.append('#PBS -m n')

        if job_tmpl.job_name:
            job_title = re.sub(r'[^a-zA-Z0-9_.-]+', '', job_tmpl.job_name)
            if not job_title or (job_title[0] not in string.ascii_letters + string.digits):
                job_title = f'j{job_title}'
            job_title = job_title[:15]
            lines.append(f'#PBS -N {job_title}')
        
        if job_tmpl.import_sys_environment:
            lines.append('#PBS -V')

        if job_tmpl.sched_output_path:
            lines.append(f'#PBS -o {job_tmpl.sched_output_path}')

        if job_tmpl.sched_join_files:
            # from qsub man page:
            # 'oe': Standard error and standard output are merged  into
            #       standard output
            # 'eo': Standard error and standard output are merged  into
            #       standard error
            # 'n' : Standard error and standard output are not merged (default)
            lines.append('#PBS -j oe')
            if job_tmpl.sched_error_path:
                print(
                    'sched_join_files is True, but sched_error_path is set in '
                    'PBSPro script; ignoring sched_error_path', flush=True
                )
        else:
            if job_tmpl.sched_error_path:
                lines.append(f'#PBS -e {job_tmpl.sched_error_path}')

        if job_tmpl.queue_name:
            lines.append(f'#PBS -q {job_tmpl.queue_name}')

        if job_tmpl.account:
            lines.append(f'#PBS -A {job_tmpl.account}')

        if job_tmpl.priority:
            lines.append(f'#PBS -p {job_tmpl.priority}')

        if not job_tmpl.job_resource:
            raise ValueError('Job resources (as the num_machines) are required for the PBSPro scheduler plugin')

        resource_lines = self._get_resource_lines(
            num_machines=job_tmpl.job_resource['num_machines'],
            num_mpiprocs_per_machine=job_tmpl.job_resource['num_mpiprocs_per_machine'],
            num_cores_per_machine=job_tmpl.job_resource['num_cores_per_machine'],
            max_memory_kb=job_tmpl.max_memory_kb,
            max_wallclock_seconds=job_tmpl.max_wallclock_seconds
        )

        lines += resource_lines

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

        lines.append('cd "$PBS_O_WORKDIR"')
        lines.append(empty_line)

        return '\n'.join(lines)

    def _get_submit_command(self, submit_script):
        """
        Return the string to execute to submit a given script.
        Args:
            submit_script: the path of the submit script
        """
        submit_command = f'qsub {submit_script}'

        return submit_command
