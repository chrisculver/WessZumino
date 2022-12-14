from collections import defaultdict

import sympy as sp
import numpy as np

import copy


def sympy_expr_to_list(expr):
    arg_list=[]

    for a in sp.preorder_traversal(expr):
        if(a!=expr):
            arg_list.append(a)

    return arg_list

