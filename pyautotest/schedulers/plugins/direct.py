'''
Date: 2021-02-23 17:13:31
LastEditors: jiyuyang
LastEditTime: 2021-04-23 16:47:06
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from pyautotest.schedulers.datastructure import NodeNumberJobResource
from pyautotest.schedulers.schedulers import Scheduler

class DirectJobResource(NodeNumberJobResource):
    pass

class DirectScheduler(Scheduler):
    """
    Support for the direct execution bypassing schedulers.
    """

    def _get_submit_script_header(self, job_tmpl):
        """
        Return the submit script header
        """

        lines = []
        empty_line = ''

        if job_tmpl.sched_output_path:
            lines.append(f'exec > {job_tmpl.sched_output_path}')

        if job_tmpl.sched_join_files:
            if job_tmpl.sched_error_path:
                print('sched_join_files is True, but sched_error_path is set; ignoring sched_error_path', flush=True)
        else:
            if job_tmpl.sched_error_path:
                lines.append(f'exec 2> {job_tmpl.sched_error_path}')
            else:
                lines.append('exec 2>&1')

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
            lines.append(f'ulimit -v {virtual_memory_kb}')
        if not job_tmpl.import_sys_environment:
            lines.append('env --ignore-environment \\')

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
        """
        submit_command = f'bash -e {submit_script} > /dev/null 2>&1 & echo $!'

        return submit_command