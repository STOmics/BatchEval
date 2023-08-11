#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/10/23 5:44 PM
# @Author  : zhangchao
# @File    : demo.py
# @Email   : zhangchao5@genomics.cn
import scanpy as sc

data = [
    sc.read_h5ad(t) for t in [
        "/media/Data/zhangchao/zhangchao/BatchEval/demo_data/stereo_olfactory_bulb_ann.h5ad",
        "/media/Data/zhangchao/zhangchao/BatchEval/demo_data/visium_olfactory_bulb_ann.h5ad"]]

from BatchEval import batch_qc

batch_qc(*data,
         qc_mode="raw",
         adjust_method="harmony",
         norm_log=True,
         is_scale=True,
         n_pcs=50,
         n_neighbors=100,
         batch_key="batch",
         position_key="X_umap",
         condition=None,
         use_rep="X_pca",
         count_key="total_counts",
         celltype_key=None,
         report_path="./output")
