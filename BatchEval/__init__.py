#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/10/23 5:10 PM
# @Author  : zhangchao
# @File    : __init__.py.py
# @Email   : zhangchao5@genomics.cn
from .BatchEval import batch_eval
from warnings import filterwarnings

filterwarnings("ignore")
__version__ = "2.0.0"
