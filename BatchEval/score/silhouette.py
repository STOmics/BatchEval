#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2023/9/6 15:10
# @Author  : zhangchao
# @File    : silhouette.py
# @Software: PyCharm
# @Email   : zhangchao5@genomics.cn
from sklearn.metrics import silhouette_score


def silhouette(adata, label_key, embed, metric="euclidean", scale=True):
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f"{embed} not in obsm")

    asw = silhouette_score(
        X=adata.obsm[embed], labels=adata.obs[label_key], metric=metric
    )
    if scale:
        asw = (asw + 1) / 2
    return asw
