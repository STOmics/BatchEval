#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/10/23 5:44 PM
# @Author  : zhangchao
# @File    : demo.py
# @Email   : zhangchao5@genomics.cn
import os
import scanpy as sc
from BatchEval import batch_eval

data = [
    sc.read_h5ad(os.path.join("./demo_data", t)) for t in sorted(os.listdir("./demo_data")) if t.endswith("h5ad")]

batch_eval(*data,
           qc_mode="raw",
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
