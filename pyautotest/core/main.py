'''
Date: 2021-03-18 14:37:53
LastEditors: jiyuyang
LastEditTime: 2021-04-25 10:35:11
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''
#TODO: refactor code with logging, warnings etc.
from pyautotest.core.convert import Convert
from pyautotest.core.run import Run
from pyautotest.core.show import Show

import argparse

def main():
    parser = argparse.ArgumentParser(prog='autotest', description='Auto-test for ABACUS')
    subparsers = parser.add_subparsers(help='sub-command help')

    # Run
    parser_run = subparsers.add_parser('run', help='way to run auto-test')
    test_group = parser_run.add_argument_group(title='Test Parameters')
    test_group.add_argument('-s', '--single', dest='single', type=str, default=None, help='local single test for one example.')
    test_group.add_argument('-b', '--batch', dest='batch', type=str, default=None, help='local batch test for many examples.' )
    
    other_group = parser_run.add_argument_group(title='Optional Parameters')
    other_group.add_argument('-l', '--local', dest='local', type=bool, default=False, help='whether to run locally. Default: False')
    other_group.add_argument('-p', '--parallel', dest='parallel', type=bool, default=False, help='whether to test in parallel. Only valid when running batch test non-locally. Default: False')
    parser_run.set_defaults(func=Run().run_cmdline)

    # Show
    parser_show = subparsers.add_parser('show', help='show information of auto-test')
    parser_show.add_argument('-l', '--lib', dest='lib', type=str, default=None, help='show example library information.')
    parser_show.add_argument('-b', '--band', dest='band', type=str, default=None, help='plot band structure and show band information.')
    parser_show.add_argument('-d', '--dos', dest='dos', type=str, default=None, help='plot density of state(DOS).')
    parser_show.set_defaults(func=Show().show_cmdline)

    # Convert
    parser_convert = subparsers.add_parser('convert', help='convert file from one format to another')
    parser_convert.add_argument('-file', '--file', dest='file', type=str, default=None, help="convert structure file format to another.")
    parser_convert.set_defaults(func=Convert().convert_cmdline)
    
    args = parser.parse_args()
    args.func(args)

#TODO: add command about showing job information, kill job
    #parser_kill = subparsers.add_parser('stop', help='stop auto-test')