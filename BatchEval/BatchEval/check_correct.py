#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2023/9/8 10:07
# @Author  : zhangchao
# @File    : check_correct.py
# @Software: PyCharm
# @Email   : zhangchao5@genomics.cn
import scanpy as sc
import pandas as pd
from collections import defaultdict
from anndata import AnnData
from typing import Union

from BatchEval.module.trainer import domain_variance_score
from BatchEval.score import get_kbet, get_lisi, silhouette
from BatchEval.test import sample_heatmap, umap_plot, joint_plot


def check_correct(merge_data: AnnData,
                  batch_key: str = "batch",
                  use_rep: str = "X_pca",
                  position_key: str = "X_umap",
                  n_neighbors: int = 100,
                  celltype_key: Union[str, None] = None,
                  report_path: str = "./",
                  gpu: Union[str, int] = "0"):
    """BatchEval Test Corrected Dataset Pipeline

    Parameters
    ----------
    merge_data: `AnnData`
        Data matrix with rows for cells and columns for genes.
    batch_key: `str`
        Label the data batches.
    use_rep: `str`
        corrected feature
    position_key: `str`
        Compute the coordinate space of the nearest neighbor.
    n_neighbors: `int`
        Calculate the nearest neighbors of a local area. default, 100.
    celltype_key: `str`
    report_path: `str`
    gpu: `str`

    Returns
    -------

    """
    data_dict = defaultdict()
    img_dict = defaultdict()

    n_batch = merge_data.obs[batch_key].cat.categories.size

    # biological evaluation
    assert use_rep in merge_data.obsm_keys()
    n_dims = merge_data.obsm[use_rep].shape[1]
    domain_df = domain_variance_score(merge_data, input_dims=n_dims, n_batch=n_batch, use_rep=use_rep,
                                      batch_key=batch_key, batch_size=4096, gpu=gpu, save_path=report_path)
    data_dict["domain"] = domain_df

    if celltype_key is None:
        sc.tl.leiden(merge_data)
        celltype_key = "leiden"

    # KBET
    reject_score, stat_mean, pvalue_mean, accept_rate = get_kbet(
        merge_data, key=batch_key, use_rep=position_key, alpha=0.05, n_neighbors=n_neighbors)
    kbet_df = pd.DataFrame(data={"Chi Mean": stat_mean, "95% P Value": pvalue_mean, "Accept Rate": accept_rate,
                                 "Reject Rate": reject_score}, index=["Score"])
    kbet_df["describe_note"] = f"Local area Chi2 test. Local sample richness test. " \
                               f"By default, UMAP coordinates are used to calculate region neighbors. " \
                               f"default neighbors: {n_neighbors}"
    # LISI
    ilisi = get_lisi(merge_data, key=batch_key, use_rep=position_key, n_neighbors=n_neighbors)
    clisi = get_lisi(merge_data, key=celltype_key, use_rep=position_key, n_neighbors=n_neighbors)
    f1_lisi = (2 * (1 - clisi) * ilisi) / (ilisi + (1 - clisi))
    lisi_df = pd.DataFrame(data={"iLISI": ilisi, f"cLISI({celltype_key})": clisi, "LISI(f1 score)": f1_lisi},
                           index=["Score"])

    # SS
    ss_batch = silhouette(merge_data, label_key=batch_key, embed=position_key, metric="euclidean", scale=True)
    ss_celltype = silhouette(merge_data, label_key=celltype_key, embed=position_key, metric="euclidean", scale=True)
    f1_ss = (2 * (1 - ss_batch) * ss_celltype) / (ss_celltype + (1 - ss_batch))

    ss_df = pd.DataFrame(data={"SS_batch": ss_batch, f"SS_{celltype_key}": ss_celltype, "SS_f1": f1_ss},
                         index=["Score"])
    data_dict["kbet"] = kbet_df
    data_dict["lisi"] = lisi_df
    data_dict["silhouette"] = ss_df

    # heatmap_gene_src = sample_heatmap(merge_data, feat_key=use_rep, metric="correlation", batch_key=batch_key)
    umap_batch_src = umap_plot(merge_data, visualize_key=batch_key)
    umap_type_src = umap_plot(merge_data, visualize_key=celltype_key)
    joint_src = joint_plot(merge_data, batch_key=batch_key, use_rep=use_rep)
    # img_dict["heatmap"] = heatmap_gene_src
    img_dict["umap_batch"] = umap_batch_src
    img_dict["umap_type"] = umap_type_src
    img_dict["joint"] = joint_src

    return data_dict, img_dict
