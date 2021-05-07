'''
Date: 2021-03-18 14:37:53
LastEditors: jiyuyang
LastEditTime: 2021-04-25 10:35:11
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''
#TODO: refactor code with typing, logging, warnings etc.
from pyautotest.core.run import Run
from pyautotest.core.show import Show

import argparse

def main():
    parser = argparse.ArgumentParser(prog='autotest', description='Auto-test for ABACUS')
    subparsers = parser.add_subparsers(help='sub-command help')

    # Run
    parser_run = subparsers.add_parser('run', help='way to run auto-test')
    test_group = parser_run.add_argument_group(title='Test Parameters')
    test_group.add_argument('-s', '--single', dest='single', type=str, default=None, help='local single test for one example, its value should be absolute path of `input.json`')
    test_group.add_argument('-b', '--batch', dest='batch', type=str, default=None, help='local batch test for many examples, its value should be absolute path of `input.json`')
    
    other_group = parser_run.add_argument_group(title='Optional Parameters')
    other_group.add_argument('-l', '--local', dest='local', type=bool, default=False, help='whether to run locally. Default: False')
    other_group.add_argument('-p', '--parallel', dest='parallel', type=bool, default=False, help='whether to test in parallel. Only valid when running batch test non-locally. Default: False')
    parser_run.set_defaults(func=Run().run_cmdline)

    # Show
    parser_show = subparsers.add_parser('show', help='show information of auto-test')
    parser_show.add_argument('-l', '--lib', dest='lib', type=str, default=None, help='show example library information, its value should be absolute path of `input.json`')
    parser_show.set_defaults(func=Show().show_cmdline)
    
    args = parser.parse_args()
    args.func(args)

#TODO: add command about showing job information, kill job
    #parser_kill = subparsers.add_parser('stop', help='stop auto-test')