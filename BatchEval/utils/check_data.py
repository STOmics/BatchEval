#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/10/23 5:14 PM
# @Author  : zhangchao
# @File    : check_data.py
# @Email   : zhangchao5@genomics.cn
from anndata import AnnData


def check_data(data: AnnData):
    # check that the count matrix contains valid inputs. More precisely, check that inputs are non-negative integers.
    if (data.X.data % 1 != 0).any().any():
        raise ValueError("The count matrix should only contain integers.")
    if (data.X.data < 0).any().any():
        raise ValueError("The count matrix should only contain non-negative values.")
