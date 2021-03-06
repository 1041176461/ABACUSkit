'''
Date: 2021-05-06 14:53:09
LastEditors: jiyuyang
LastEditTime: 2021-05-06 15:01:31
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''
import typing
from os import PathLike

from abacuskit.schedulers.data import Code

Dict_str_str = typing.Dict[str, str]
Dict_str_list = typing.Dict[str, list]
Dict_str_int = typing.Dict[str, int]
Dict_str_float = typing.Dict[str, float]
Dict_Tuple_Dict = typing.Dict[typing.Tuple[str, str], dict]
str_PathLike = typing.Union[str, PathLike]
muti_Path = typing.Sequence[str_PathLike]
muti_Dict = typing.Sequence[typing.Dict]
Command = typing.Union[str, Code]
List_Command = typing.Sequence[Command]
Return_2 = typing.Tuple[list, list]
Int_List = typing.Union[int, list]
