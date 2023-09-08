#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2023/9/6 11:15
# @Author  : zhangchao
# @File    : harmony.py
# @Software: PyCharm
# @Email   : zhangchao5@genomics.cn
import scanpy as sc
from anndata import AnnData

from BatchEval.utils import pca_lowrank


def harmony_integration(merge_data: AnnData, batch_key: str = "batch", n_pcs: int = 50, use_rep: str = "X_pca"):
    if use_rep not in merge_data.obsm_keys():
        pca_lowrank(merge_data, n_component=n_pcs)
        use_rep = "X_pca"

    sc.external.pp.harmony_integrate(merge_data, key=batch_key, basis=use_rep)
