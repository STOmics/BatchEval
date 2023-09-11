#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2023/9/7 16:27
# @Author  : zhangchao
# @File    : check_raw.py
# @Software: PyCharm
# @Email   : zhangchao5@genomics.cn
import scanpy as sc
import pandas as pd
from collections import defaultdict
from anndata import AnnData
from typing import Union

from BatchEval.module.trainer import domain_variance_score
from BatchEval.score import get_kbet, get_lisi
from BatchEval.score.silhouette import silhouette
from BatchEval.test import variance_test, qq_plot, description_data, cdf_plot, kernel_plot, var_mean_plot, \
    sample_heatmap, umap_plot, joint_plot
from BatchEval.utils import check_data, pca_lowrank


def check_raw(
        *data: AnnData,
        norm_log: bool = False,
        is_scale: bool = False,
        n_pcs: int = 50,
        n_neighbors: int = 100,
        batch_key: str = "batch",
        position_key: str = "X_umap",
        condition: Union[str, list, None] = None,
        count_key: str = "total_counts",
        celltype_key: Union[str, None] = None,
        report_path: str = "./",
        gpu: Union[str, int] = "0"):
    """BatchEval Test Raw Dataset Pipeline

    Parameters
    ----------
    data: `Anndata`
        Data matrix with rows for cells and columns for genes.
    norm_log: `bool`
        Whether to preprocess data. `sc.pp.normalization()`, `sc.pp.log1p()`
    is_scale: `bool`
        Whether to preprocess data. `sc.pp.scale()`
    n_pcs: `int`
        Number of principal components retained in PCA. default, 50.
    n_neighbors: `int`
        Calculate the nearest neighbors of a local area. default, 100.
    batch_key: `str`
        Label the data batches.
    position_key: `str`
        Compute the coordinate space of the nearest neighbor.
    condition: `str`
        Label the experimental conditions. By default, the experimental conditions for each data are different.
    count_key: `str`
    celltype_key: `str`
    report_path: `str`
    gpu: `str`

    Returns
    -------

    """
    data_dict = defaultdict()
    img_dict = defaultdict()
    merge_data = AnnData.concatenate(*data, batch_key=batch_key)
    check_data(merge_data)
    n_batch = merge_data.obs[batch_key].cat.categories.size

    if count_key not in merge_data.obs_keys():
        sc.pp.calculate_qc_metrics(merge_data, inplace=True)
    if norm_log:
        sc.pp.normalize_total(merge_data, target_sum=1e4)
        sc.pp.log1p(merge_data)
    if is_scale:
        sc.pp.scale(merge_data, zero_center=False, max_value=10)
    pca_lowrank(merge_data, use_rep=None, n_component=n_pcs)

    sc.pp.neighbors(merge_data)
    sc.tl.umap(merge_data)

    # statistical evaluation
    F_test_df, boxplot_src = variance_test(merge_data, batch_key=batch_key, test_key=count_key)
    data_dict["f-test"] = F_test_df
    img_dict["box"] = boxplot_src

    qq_srcs, ks_static_df = qq_plot(merge_data, batch_key=batch_key, test_key=count_key)
    data_dict["ks-test"] = ks_static_df
    img_dict["qq"] = qq_srcs

    describe_df, confound_df = description_data(merge_data, condition=condition)
    data_dict["describe"] = describe_df
    data_dict["confound"] = confound_df

    cdf_src = cdf_plot(merge_data, batch_key=batch_key, use_key=count_key)
    kernel_src = kernel_plot(merge_data, batch_key=batch_key, test_key=count_key)
    vm_src = var_mean_plot(merge_data, batch_key=batch_key)
    img_dict["cdf"] = cdf_src
    img_dict["kernel"] = kernel_src
    img_dict["var-mean"] = vm_src

    # biological evaluation
    domain_df = domain_variance_score(merge_data, input_dims=n_pcs, n_batch=n_batch, use_rep="X_pca",
                                      batch_key=batch_key, batch_size=4096, gpu=gpu, save_path=report_path)
    data_dict["domain"] = domain_df

    if celltype_key is None:
        sc.tl.leiden(merge_data)
        celltype_key = "leiden"

    assert celltype_key in merge_data.obs_keys()

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
    lisi_df = pd.DataFrame(data={"iLISI": ilisi, f"cLISI({celltype_key})": clisi, "LISI(f1 score)": f1_lisi}, index=["Score"])

    # SS
    ss_batch = silhouette(merge_data, label_key=batch_key, embed=position_key, metric="euclidean", scale=True)
    ss_celltype = silhouette(merge_data, label_key=celltype_key, embed=position_key, metric="euclidean", scale=True)
    f1_ss = (2 * (1 - ss_batch) * ss_celltype) / (ss_celltype + (1 - ss_batch))

    ss_df = pd.DataFrame(data={"SS_batch": ss_batch, f"SS_{celltype_key}": ss_celltype, "SS_f1": f1_ss}, index=["Score"])
    data_dict["kbet"] = kbet_df
    data_dict["lisi"] = lisi_df
    data_dict["silhouette"] = ss_df

    # heatmap_gene_src = sample_heatmap(merge_data, feat_key="X_pca", metric="correlation", batch_key=batch_key)
    umap_batch_src = umap_plot(merge_data, visualize_key=batch_key)
    umap_type_src = umap_plot(merge_data, visualize_key=celltype_key)
    joint_src = joint_plot(merge_data, batch_key=batch_key, use_rep="X_pca")
    # img_dict["heatmap"] = heatmap_gene_src
    img_dict["umap_batch"] = umap_batch_src
    img_dict["umap_type"] = umap_type_src
    img_dict["joint"] = joint_src

    return merge_data, data_dict, img_dict
