'''
Date: 2021-02-11 17:35:09
LastEditors: jiyuyang
LastEditTime: 2021-04-28 14:39:14
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

import abc

from abacuskit.schedulers.data import CodeRunMode, JobDefaultFields


class Scheduler(abc.ABC):

    def __str__(self):
        return self.__class__.__name__

    def get_submit_script(self, job_tmpl):
        """Return the submit script as a string"""

        assert isinstance(job_tmpl, JobDefaultFields)

        empty_line = ''
        script_lines = []

        if job_tmpl.shebang:
            script_lines.append(job_tmpl.shebang)
        elif job_tmpl.shebang == '':
            script_lines.append(job_tmpl.shebang)
        elif job_tmpl.shebang is None:
            script_lines.append('#!/bin/bash')
        else:
            raise ValueError(f'Invalid shebang set: {job_tmpl.shebang}')
        script_lines.append(self._get_submit_script_header(job_tmpl))
        script_lines.append(empty_line)

        if job_tmpl.prepend_text:
            script_lines.append(job_tmpl.prepend_text)
            script_lines.append(empty_line)

        script_lines.append(self._get_run_line(
            job_tmpl.codes_info, job_tmpl.codes_run_mode))
        script_lines.append(empty_line)

        if job_tmpl.append_text:
            script_lines.append(job_tmpl.append_text)
            script_lines.append(empty_line)

        footer = self._get_submit_script_footer(job_tmpl)
        if footer:
            script_lines.append(footer)
            script_lines.append(empty_line)

        return '\n'.join(script_lines)

    @abc.abstractmethod
    def _get_submit_script_header(self, job_tmpl):
        """Return the submit script header"""

    def _get_submit_script_footer(self, job_tmpl):
        """Return the submit script final part"""
        return None

    def _get_run_line(self, codes_info, codes_run_mode):
        """Return a string with the line to execute a specific code with specific arguments.

        :parameter codes_info: a list of `Code` objects. Each contains the information needed to run the code.
        :parameter codes_run_mode: instance of `CodeRunMode` contains the information on how to launch the multiple codes.
        :return: string with format: [executable] [args] {[ < stdin ]} {[ < stdout ]} {[2>&1 | 2> stderr]}
        """

        list_of_runlines = [code.run_line() for code in codes_info]

        if codes_run_mode == CodeRunMode.PARALLEL:
            list_of_runlines.append('wait\n')
            return ' &\n\n'.join(list_of_runlines)

        if codes_run_mode == CodeRunMode.SERIAL:
            return '\n\n'.join(list_of_runlines)

        raise NotImplementedError('Unrecognized code run mode')

    @abc.abstractmethod
    def _get_submit_command(self, submit_script):
        """Return the string to execute to submit a given script.

        :param submit_script: the path of the submit script relative to the working directory.
        :return: the string to execute to submit a given script.
        """
