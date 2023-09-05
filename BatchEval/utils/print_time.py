#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/10/23 5:17 PM
# @Author  : zhangchao
# @File    : print_time.py
# @Email   : zhangchao5@genomics.cn
import time


def print_time(function):
    def func_time(*args, **kwargs):
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} BatchEval Pipeline Starting")
        t0 = time.time()
        res = function(*args, **kwargs)
        t1 = time.time()
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} BatchEval Pipeline Done!")
        t = t1 - t0
        print(f"Total Running Time: {t // 60:}min {t % 60:.4f}s")
        return res

    return func_time
