#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2023/9/6 11:16
# @Author  : zhangchao
# @File    : bbknn.py
# @Software: PyCharm
# @Email   : zhangchao5@genomics.cn
import scanpy as sc
from anndata import AnnData

from BatchEval.utils import pca_lowrank


def bbknn_integration(merge_data: AnnData, batch_key: str = "batch", n_pcs: int = 50, use_rep: str = "X_pca"):
    if use_rep not in merge_data.obsm_keys():
        pca_lowrank(merge_data, n_component=n_pcs)
        use_rep = "X_pca"
    sc.external.pp.bbknn(merge_data, batch_key=batch_key, use_rep=use_rep)
