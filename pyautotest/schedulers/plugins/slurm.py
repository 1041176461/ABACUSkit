'''
Date: 2021-02-19 11:16:01
LastEditors: jiyuyang
LastEditTime: 2021-04-23 16:47:57
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.schedulers.schedulers import Scheduler
from pyautotest.schedulers.datastructure import NodeNumberJobResource


class SlurmJobResource(NodeNumberJobResource):
    """Class for SLURM job resource"""

    @classmethod
    def validate_resources(cls, **kwargs):
        """Validate the resources against the job resource class of this scheduler.
        
        This extends the base class validator to check that the `num_cores_per_machine` are a multiple of
        `num_cores_per_mpiproc` and/or `num_mpiprocs_per_machine`.

        :param kwargs: dictionary of values to define the job resources
        :return: defaultdict with the parsed parameters populated
        :raises ValueError: if the resources are invalid or incomplete
        """
        resources = super().validate_resources(**kwargs)

        # In this plugin we never used num_cores_per_machine so if it is not defined it is OK.
        if resources["num_cores_per_machine"] is not None and resources["num_cores_per_mpiproc"] is not None:
            if resources["num_cores_per_machine"] != resources["num_cores_per_mpiproc"] * resources["num_mpiprocs_per_machine"]:
                raise ValueError(
                    '`num_cores_per_machine` must be equal to `num_cores_per_mpiproc * num_mpiprocs_per_machine` and in'
                    ' particular it should be a multiple of `num_cores_per_mpiproc` and/or `num_mpiprocs_per_machine`'
                )

        elif resources["num_cores_per_machine"] is not None:
            if resources["num_cores_per_machine"] < 1:
                raise ValueError('num_cores_per_machine must be greater than or equal to one.')

            resources["num_cores_per_mpiproc"] = (resources["num_cores_per_machine"] / resources["num_mpiprocs_per_machine"])
            if int(resources["num_cores_per_mpiproc"]) != resources["num_cores_per_mpiproc"]:
                raise ValueError(
                    '`num_cores_per_machine` must be equal to `num_cores_per_mpiproc * num_mpiprocs_per_machine` and in'
                    ' particular it should be a multiple of `num_cores_per_mpiproc` and/or `num_mpiprocs_per_machine`'
                )
            resources["num_cores_per_mpiproc"] = int(resources["num_cores_per_mpiproc"])

        return resources

class SlurmScheduler(Scheduler):
    """
    Support for the SLURM scheduler (http://slurm.schedmd.com/).
    """

    def _get_submit_script_header(self, job_tmpl):
        """
        Return the submit script header, using the parameters from the job_tmpl.
        """
        import re
        import string

        empty_line = ''

        lines=[]
        if job_tmpl.submit_as_hold:
            lines.append('#SBATCH -H')

        if job_tmpl.rerunnable:
            lines.append('#SBATCH --requeue')
        else:
            lines.append('#SBATCH --no-requeue')

        if job_tmpl.email:
            lines.append(f'#SBATCH --mail-user={job_tmpl.email}')

        if job_tmpl.email_on_started:
            lines.append('#SBATCH --mail-type=BEGIN')
        if job_tmpl.email_on_terminated:
            lines.append('#SBATCH --mail-type=FAIL')
            lines.append('#SBATCH --mail-type=END')

        if job_tmpl.job_name:
            job_title = re.sub(r'[^a-zA-Z0-9_.-]+', '', job_tmpl.job_name)
            if not job_title or (job_title[0] not in string.ascii_letters + string.digits):
                job_title = f'j{job_title}'
            job_title = job_title[:128]
            lines.append(f'#SBATCH --job-name="{job_title}"')

        if job_tmpl.import_sys_environment:
            lines.append('#SBATCH --get-user-env')

        if job_tmpl.sched_output_path:
            lines.append(f'#SBATCH --output={job_tmpl.sched_output_path}')

        if job_tmpl.sched_join_files:
            # By  default both standard output and standard error are directed
            # to a file of the name "slurm-%j.out", where the "%j" is replaced
            # with  the  job  allocation  number.
            # See that this automatic redirection works also if
            # I specify a different --output file
            if job_tmpl.sched_error_path:
                print(
                    'sched_join_files is True, but sched_error_path is set in '
                    'SLURM script; ignoring sched_error_path', flush=True
                )
        else:
            if job_tmpl.sched_error_path:
                lines.append(f'#SBATCH --error={job_tmpl.sched_error_path}')
            else:
                # To avoid automatic join of files
                lines.append('#SBATCH --error=slurm-%j.err')

        if job_tmpl.queue_name:
            lines.append(f'#SBATCH --partition={job_tmpl.queue_name}')

        if job_tmpl.account:
            lines.append(f'#SBATCH --account={job_tmpl.account}')

        if job_tmpl.qos:
            lines.append(f'#SBATCH --qos={job_tmpl.qos}')

        if job_tmpl.priority:
            #  Run the job with an adjusted scheduling priority  within  SLURM.
            #  With no adjustment value the scheduling priority is decreased by
            #  100. The adjustment range is from -10000 (highest  priority)  to
            #  10000  (lowest  priority).
            lines.append(f'#SBATCH --nice={job_tmpl.priority}')

        if not job_tmpl.job_resource:
            raise ValueError('Job resources (as the num_machines) are required for the SLURM scheduler plugin')

        lines.append(f'#SBATCH --nodes={job_tmpl.job_resource["num_machines"]}')
        if job_tmpl.job_resource["num_mpiprocs_per_machine"]:
            lines.append(f'#SBATCH --ntasks-per-node={job_tmpl.job_resource["num_mpiprocs_per_machine"]}')

        if job_tmpl.job_resource["num_cores_per_mpiproc"]:
            lines.append(f'#SBATCH --cpus-per-task={job_tmpl.job_resource["num_cores_per_mpiproc"]}')

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
            days = tot_secs // 86400
            tot_hours = tot_secs % 86400
            hours = tot_hours // 3600
            tot_minutes = tot_hours % 3600
            minutes = tot_minutes // 60
            seconds = tot_minutes % 60
            if days == 0:
                lines.append(f'#SBATCH --time={hours:02d}:{minutes:02d}:{seconds:02d}')
            else:
                lines.append(f'#SBATCH --time={days:d}-{hours:02d}:{minutes:02d}:{seconds:02d}')

        # It is the memory per node, not per cpu!
        if job_tmpl.max_memory_kb:
            try:
                virtual_memory_kb = int(job_tmpl.max_memory_kb)
                if virtual_memory_kb <= 0:
                    raise ValueError
            except ValueError:
                raise ValueError(
                    'max_memory_kb must be '
                    "a positive integer (in kB)! It is instead '{}'"
                    ''.format((job_tmpl.MaxMemoryKb))
                )
            # --mem: Specify the real memory required per node in MegaBytes.
            lines.append(f'#SBATCH --mem={virtual_memory_kb // 1024}')

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

        lines.append(empty_line)

        return '\n'.join(lines)

    def _get_submit_command(self, submit_script):
        """
        Return the string to execute to submit a given script.

        Args:
            submit_script: the path of the submit script
        """
        submit_command = f'sbatch {submit_script}'

        return submit_command
