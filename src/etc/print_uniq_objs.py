#!/usr/bin/env python
"""
print_uniq_objs.py

usage:
    print_uniq_objs.py sourcefiles...
"""
import os
import sys

OBJ_SUFFIX = '.o'


def trunkname(name):
    """
    path1/path2/aaa.bbb.ccc.F90 => aaa
    """
    return os.path.basename(name).split('.', 1)[0]


def uniq_with_trunkname(original_names):
    names = {}
    for x in original_names:
        names[trunkname(x)] = x

    return list(names.values())


if __name__ == '__main__':
    for x in uniq_with_trunkname(sys.argv[1:]):
        print(os.path.splitext(x)[0] + OBJ_SUFFIX)
