#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/10/23 5:38 PM
# @Author  : zhangchao
# @File    : check_raw.py
# @Email   : zhangchao5@genomics.cn
import os
import os.path as osp
import pkgutil
import pwd
import time
import scanpy as sc
import pandas as pd
from typing import Union
from anndata import AnnData
from lxml import etree
from collections import defaultdict

from BatchEval.module.trainer import domain_variance_score
from BatchEval.test import description_data, variance_test, cdf_plot, kernel_plot, var_mean_plot, qq_plot, \
    distribution_fitting, metric_score, sample_heatmap, umap_plot, joint_plot
from BatchEval.utils import check_data, pca_lowrank, embed_text, embed_tabel, embed_table_imgs


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
        gpu: Union[str, int] = "0") -> dict:
    """BatchEval Test Raw Dataset Pipeline

    Parameters
    -----------------
    *data: 'Anndata'
        Data matrix with rows for cells and columns for genes.
    norm_log: 'bool'
        Whether to preprocess data. 'sc.pp.normalization()', 'sc.pp.log1p()'
    is_scale: 'bool'
         Whether to preprocess data. 'sc.pp.scale()'
    n_pcs: 'int'
        Number of principal components retained in PCA. default, 50.
    n_neighbors: 'int'
        Calculate the nearest neighbors of a local area. default, 100.
    batch_key: 'str'
        Label the data batches.
    position_key: 'str'
        Compute the coordinate space of the nearest neighbor.
    condition: 'str, list, None'
        Label the experimental conditions. By default, the experimental conditions for each data are different.
    count_key: 'str'
    celltype_key: 'str'
    report_path: 'str'
    gpu: 'str', 'int'

    Return
    -----------------
    output_dict: 'dict'
    """
    output_dict = defaultdict()
    merge_data = AnnData.concatenate(*data, batch_key=batch_key)
    n_batch = merge_data.obs[batch_key].cat.categories.size

    check_data(merge_data)

    if merge_data.raw is None:
        merge_data.raw = merge_data
    else:
        merge_data = merge_data.raw.to_adata()
    if count_key not in merge_data.obs_keys():
        sc.pp.calculate_qc_metrics(merge_data, inplace=True)
    if norm_log:
        sc.pp.normalize_total(merge_data, target_sum=1e4)
        sc.pp.log1p(merge_data)
    if is_scale:
        sc.pp.scale(merge_data, zero_center=False, max_value=10)
    pca_lowrank(merge_data, use_rep=None, n_component=n_pcs)
    output_dict["dataset"] = merge_data.copy()

    sc.pp.neighbors(merge_data)
    sc.tl.umap(merge_data)

    domain_df = domain_variance_score(merge_data,
                                      input_dims=n_pcs,
                                      n_batch=n_batch,
                                      use_rep="X_pca",
                                      batch_key=batch_key,
                                      batch_size=4096,
                                      gpu=gpu,
                                      save_path=report_path)
    output_dict["table"] = {"domain": domain_df}

    describe_df, confound_df = description_data(merge_data, condition=condition)
    output_dict["table"].update({"describe": describe_df, "confound": confound_df})

    F_test_df, boxplot_src = variance_test(merge_data, batch_key=batch_key, test_key=count_key)
    output_dict["table"]["f-test"] = F_test_df
    output_dict["imgs"] = {"box": boxplot_src}

    cdf_src = cdf_plot(merge_data, batch_key=batch_key, use_key=count_key)
    kernel_src = kernel_plot(merge_data, batch_key=batch_key, test_key=count_key)
    vm_src = var_mean_plot(merge_data, batch_key=batch_key)
    output_dict["imgs"]["cdf"] = cdf_src
    output_dict["imgs"]["kernel"] = kernel_src
    output_dict["imgs"]["var-mean"] = vm_src

    qq_srcs, ks_static_df = qq_plot(merge_data, batch_key=batch_key, test_key=count_key)
    output_dict["table"]["ks"] = ks_static_df
    output_dict["imgs"]["qq"] = qq_srcs
    dist_srcs = distribution_fitting(merge_data, batch_key=batch_key, fit_key=count_key)
    output_dict["imgs"]["dist"] = dist_srcs

    metric_dict = metric_score(merge_data, n_neighbor=n_neighbors, batch_key=batch_key, metric_pos=position_key,
                               celltype_key=celltype_key)
    output_dict["table"].update(metric_dict)

    heatmap_gene_src = sample_heatmap(merge_data, feat_key="X_pca", metric="correlation", batch_key=batch_key)
    umap_batch_src = umap_plot(merge_data, visualize_key=batch_key)
    if celltype_key is not None:
        umap_type_src = umap_plot(merge_data, visualize_key=celltype_key)
        output_dict["imgs"]["umap_type"] = umap_type_src
    joint_srcs = joint_plot(merge_data, batch_key=batch_key, use_rep="X_pca")
    output_dict["imgs"]["heatmap"] = heatmap_gene_src
    output_dict["imgs"]["umap_batch"] = umap_batch_src
    output_dict["imgs"]["joint"] = joint_srcs

    batch_score = output_dict["table"]["domain"].iloc[0]["Accept Rate"] * \
                  output_dict["table"]["confound"].iloc[0]["Cramer's V Coefficient"] * \
                  output_dict["table"]["kbet_df"].iloc[0]["Accept Rate"] * \
                  output_dict["table"]["lisi_df"].iloc[0]["LISI Mean"] / n_batch

    is_scale = False
    if output_dict["table"]["kbet_df"].iloc[0]["95% P Value"] < 0.05:
        threshold = 0.05
    elif output_dict["table"]["kbet_df"].iloc[0]["95% P Value"] < 0.1:
        is_scale = True
        threshold = 0.1
    else:
        threshold = 0.1

    if batch_score < threshold:
        conclusion = "Need to do batch effect removal."
    elif is_scale:
        conclusion = "The data difference is small. Recommend approach: 'z-score'."
    else:
        conclusion = "Don't have to do batch effects removal."

    summary_df = pd.DataFrame(
        data={"score": f"{batch_score:.4f}", "adaptive threshold": threshold, "conclusion": conclusion},
        index=["result"])

    output_dict["table"]["summary"] = summary_df
    generate_report(data_dict=output_dict, save_path=report_path)

    return output_dict


def generate_report(data_dict: dict,
                    save_path: str) -> None:
    """Generate BatchEval Report

    Parameters
    -----------------
    data_dict: 'dict'
    save_path: 'str'
    """
    html = etree.HTML(pkgutil.get_data("BatchEval", "report_template_raw.html").decode())

    # -------- set username & run time --------
    embed_text(html, pos="h4", name="username", text=f"Report By: {pwd.getpwuid(os.getuid())[0]}")
    embed_text(html, pos="h5", name="runtime",
               text=f"Report Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")

    # -------- insert table --------
    embed_tabel(data_dict["table"]["describe"], html, pos="h4", name="describe", is_round=False)
    embed_tabel(data_dict["table"]["f-test"], html, pos="h4", name="f-test")
    embed_tabel(data_dict["table"]["ks"], html, pos="h4", name="ks-test")
    embed_tabel(data_dict["table"]["domain"], html, pos="h4", name="domain-score")
    embed_tabel(data_dict["table"]["confound"], html, pos="h4", name="confound")

    embed_tabel(data_dict["table"]['kbet_df'], html, pos="h4", name="kbet")
    embed_tabel(data_dict["table"]['lisi_df'], html, pos="h4", name="lisi")
    embed_tabel(data_dict["table"]['ksim_df'], html, pos="h4", name="ksim")
    embed_tabel(data_dict["table"]['summary'], html, pos="h4", name="conclusion")

    # -------- insert images --------
    src_dict = {
        "Kernel Distribution Curve of UMICount Total": data_dict["imgs"]["kernel"],
        "CDF Curve of UMICount Total": data_dict["imgs"]["cdf"]
    }
    embed_table_imgs(src_dict, tree=html, pos="div", class_name="curve")

    src_dict = {
        "Box-plot of UMICount Total": data_dict["imgs"]["box"],
        "UMICount Scatter Plot for Each Spot": data_dict["imgs"]["var-mean"]
    }
    embed_table_imgs(src_dict, tree=html, pos="div", class_name="plot")

    src_dict = {
        "HeatMap": data_dict["imgs"]["heatmap"],
        "Joint": data_dict["imgs"]["joint"],
        "UMAP-Batch": data_dict["imgs"]["umap_batch"]
    }
    if "umap_type" in data_dict["imgs"].keys():
        src_dict.update({"UMAP-Type": data_dict["imgs"]["umap_type"]})
    embed_table_imgs(buffer_dict=src_dict, tree=html, pos="h4", class_name="sample")
    embed_table_imgs(buffer_dict=data_dict["imgs"]["qq"], tree=html, pos="h4", class_name="dist-qq")
    embed_table_imgs(buffer_dict=data_dict["imgs"]["dist"], tree=html, pos="h4", class_name="dist-norm")

    tree = etree.ElementTree(html)
    os.makedirs(save_path, exist_ok=True)
    tree.write(osp.join(save_path, "BatchEval_report_raw.html"))
    print(
        f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} The 'BatchEval_report_raw.html' has been saved to {save_path}")
