#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/10/23 5:44 PM
# @Author  : zhangchao
# @File    : demo.py
# @Email   : zhangchao5@genomics.cn
import os
import scanpy as sc
from anndata import AnnData
from BatchEval import batch_eval

data = [
    sc.read_h5ad(os.path.join("./demo_data", t)) for t in sorted(os.listdir("./demo_data")) if t.endswith("h5ad")]

merge_data = AnnData.concatenate(*data)
sc.pp.normalize_total(merge_data)
sc.pp.log1p(merge_data)
sc.pp.combat(merge_data)

batch_eval(*data,
           norm_log=True,
           is_scale=False,
           n_pcs=50,
           n_neighbors=15,
           batch_key="batch",
           position_key="X_umap",
           condition=None,
           count_key="total_counts",
           celltype_key=None,
           external_list=[{
               "correct_data": merge_data,
               "re_pca": True,
               "re_neigh": True,
               "use_rep": "X_pca",
               "re_umap": True,
               "batch_key": "batch"},],
           report_path="./output/mouse_ob")
