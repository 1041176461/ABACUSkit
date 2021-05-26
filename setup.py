'''
Date: 2021-03-17 21:12:53
LastEditors: jiyuyang
LastEditTime: 2021-04-27 17:34:16
Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
'''

from setuptools import setup, find_packages

if __name__ == "__main__":
    setup(
        name = 'abacuskit',
        version = '1.0.0',
        packages = find_packages(),
        description = 'Autotest and Postprocessing tool for ABACUS code',
        author = 'jiyuyang',
        author_email = 'jiyuyang@mail.ustc.edu.cn',
        url = 'None',
        entry_points={'console_scripts': ['abacuskit=abacuskit.core.main:main']}
    )
