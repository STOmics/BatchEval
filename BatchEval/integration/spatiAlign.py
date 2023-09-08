#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2023/9/6 10:32
# @Author  : zhangchao
# @File    : spatiAlign.py
# @Software: PyCharm
# @Email   : zhangchao5@genomics.cn

from spatialign import Spatialign


def spatiAlign_integration(merge_data, batch_key, save_path, tau1=0.2, tau2=1., tau3=0.5):
    model = Spatialign(
        merge_data, batch_key=batch_key, is_hvg=False, is_reduce=False, n_pcs=100, n_hvg=2000, n_neigh=15,
        is_undirected=True, latent_dims=100, is_verbose=False, seed=42, gpu=0, save_path=save_path
    )
    model.train(lr=1e-3, max_epoch=500, alpha=0.5, patient=15, tau1=tau1, tau2=tau2, tau3=tau3)

    corrected_data = model.alignment()

    return corrected_data
