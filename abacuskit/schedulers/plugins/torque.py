'''
Date: 2021-02-21 00:15:49
LastEditors: jiyuyang
LastEditTime: 2021-02-21 00:24:54
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''
from .pbsbaseclasses import PbsBaseClass


class TorqueScheduler(PbsBaseClass):
    """
    Subclass to support the Torque scheduler..
    """

    def _get_resource_lines(
        self, num_machines, num_mpiprocs_per_machine, num_cores_per_machine, max_memory_kb, max_wallclock_seconds
    ):
        """
        Return the lines for machines, memory and wallclock relative
        to pbspro.
        """
        return_lines = []

        select_string = f'nodes={num_machines}'
        if num_cores_per_machine:
            select_string += f':ppn={num_cores_per_machine}'
        elif num_mpiprocs_per_machine:
            # if num_cores_per_machine is not defined then use
            # num_mpiprocs_per_machine
            select_string += f':ppn={num_mpiprocs_per_machine}'

        if max_wallclock_seconds is not None:
            try:
                tot_secs = int(max_wallclock_seconds)
                if tot_secs <= 0:
                    raise ValueError
            except ValueError:
                raise ValueError(
                    'max_wallclock_seconds must be '
                    "a positive integer (in seconds)! It is instead '{}'"
                    ''.format(max_wallclock_seconds)
                )
            hours = tot_secs // 3600
            tot_minutes = tot_secs % 3600
            minutes = tot_minutes // 60
            seconds = tot_minutes % 60
            # There is always something before, at least the total #
            # of nodes
            select_string += f',walltime={hours:02d}:{minutes:02d}:{seconds:02d}'

        if max_memory_kb:
            try:
                virtual_memory_kb = int(max_memory_kb)
                if virtual_memory_kb <= 0:
                    raise ValueError
            except ValueError:
                raise ValueError(
                    'max_memory_kb must be '
                    "a positive integer (in kB)! It is instead '{}'"
                    ''.format((max_memory_kb))
                )
            select_string += f',mem={virtual_memory_kb}kb'

        return_lines.append(f'#PBS -l {select_string}')
        return return_lines
