'''
Date: 2021-05-10 11:40:58
LastEditors: jiyuyang
LastEditTime: 2021-05-10 11:40:58
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

import cProfile, pstats
from typing import Union

def profile(column:str="time", list:Union[int, float, str]=5):
    def _profile(function):
        def _profile(*args, **kwargs):
            profiler = cProfile.Profile()
            profiler.enable()
            result = function(*args, **kwargs)
            profiler.disable()
            p = pstats.Stats(profiler)
            p.sort_stats(column).print_stats(list)
            return result
        return _profile
    return _profile