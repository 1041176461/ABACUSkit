'''
Date: 2021-02-13 18:16:34
LastEditors: jiyuyang
LastEditTime: 2021-02-13 18:38:21
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from .pbsbaseclasses import PbsBaseClass


class PbsproScheduler(PbsBaseClass):
    """
    Subclass to support the PBSPro scheduler
    (http://www.pbsworks.com/).
    """

    def _get_resource_lines(self, num_machines, num_mpiprocs_per_machine, num_cores_per_machine, max_memory_kb, max_wallclock_seconds):
        """
        Return the lines for machines, memory and wallclock relative to pbspro.
        """
        return_lines = []

        select_string = f'select={num_machines}'
        if num_mpiprocs_per_machine:
            select_string += f':mpiprocs={num_mpiprocs_per_machine}'
        if num_cores_per_machine:
            select_string += f':ppn={num_cores_per_machine}'

        if max_wallclock_seconds is not None:
            tot_secs = int(max_wallclock_seconds)
            if tot_secs <= 0:
                raise ValueError(
                    '`max_wallclock_seconds` must a positive integer!')
            hours = tot_secs // 3600
            tot_minutes = tot_secs % 3600
            minutes = tot_minutes // 60
            seconds = tot_minutes % 60
            return_lines.append(
                f'#PBS -l walltime={hours:02d}:{minutes:02d}:{seconds:02d}')

        if max_memory_kb:
            virtual_memory_kb = int(max_memory_kb)
            if virtual_memory_kb <= 0:
                raise ValueError('`max_memory_kb` must a positive integer!')
            select_string += f':mem={virtual_memory_kb}kb'

        return_lines.append(f'#PBS -l {select_string}')
        return return_lines
